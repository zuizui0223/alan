#!/usr/bin/env python3
"""Exploratory macro-analysis of recorded GloBI moth-visitor plants and historic ALAN.

Input: SOURCE-VERIFIED full93m same-year GSHHG island cell-year/WCVP accepted
recorded-species numerator AND denominator, with PRE-PHOTO 1992-2024 PCNL.
No flower ROI, no imputed biological absence, no effective nocturnal
pollination inference, no causal plant-decline claims.

Primary: within-island, year-fixed-effect linear probability association
of *recorded GloBI moth-positive accepted-taxon fraction* with harmonized
mean log(1+PCNL), adjustment for photo observation effort, each island equal.
Site-equal, richer-site adjustment and top-island deletion are sensitivities.
Never tune sampling or models to obtain a significant p-value.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

YEARS=(2019,2021,2023,2024,2025)
LIGHT="pcnl_1992_to_previous_year_mean_log1p"
REQUIRED={"island","latbin","lonbin","year",
          "all_accepted_recorded_species",
          "moth_visitor_positive_recorded_species",
          "all_accepted_observation_uuids","all_accepted_photo_ids",
          "moth_positive_fraction_OF_RECORDED_WCVP_SPECIES",
          "pcnl_years_expected","pcnl_years_observed",
          "historic_ntl_complete",LIGHT}
SCHEMA="ALAN_real_full93m_island_GloBI_positive_fraction_historic_PCNLPrior_two_way_FE_v1"


def require(condition,message):
    if not condition:raise ValueError("ALAN_PCNL_MOTH_MACRO: "+str(message))


def prepare(frame, obs_min=5, species_min=3, light_range_min=0.02):
    require(REQUIRED.issubset(frame.columns),
            "Input is not signed full93m allflora numerator/denominator + PCNL data")
    x=frame.copy()
    x["island"]=x.island.astype(str)
    for key in ("year","latbin","lonbin","all_accepted_recorded_species",
                "moth_visitor_positive_recorded_species",
                "all_accepted_observation_uuids","all_accepted_photo_ids",
                "pcnl_years_expected","pcnl_years_observed",LIGHT):
        x[key]=pd.to_numeric(x[key],errors="coerce")
    require(x[["island","latbin","lonbin","year"]].duplicated().sum()==0,
            "Double counted real photo island/site/year")
    require(set(x.year.dropna().astype(int)).issubset(set(YEARS)) and
            x.year.notna().all(),"Wrong original photo year")
    require(x.all_accepted_recorded_species.ge(
                x.moth_visitor_positive_recorded_species).all() and
            x.all_accepted_observation_uuids.le(x.all_accepted_photo_ids).all(),
            "Positive evidence or observation denominator exceeds original source")
    ratio=(x.moth_visitor_positive_recorded_species /
           x.all_accepted_recorded_species.replace(0,np.nan))
    require(np.allclose(ratio,x.moth_positive_fraction_OF_RECORDED_WCVP_SPECIES,
                        equal_nan=True,atol=1e-9,rtol=0),
            "Saved GloBI positive documented-species fraction not source-derived")
    valid=(x.all_accepted_observation_uuids.ge(obs_min) &
           x.all_accepted_recorded_species.ge(species_min) &
           x.historic_ntl_complete.astype(str).str.lower().eq("true") &
           x.pcnl_years_observed.eq(x.pcnl_years_expected) &
           x.pcnl_years_expected.eq(x.year-1992) &
           np.isfinite(x[LIGHT].to_numpy(float)) &
           x[LIGHT].ge(0))
    x=x.loc[valid].copy()
    x["frac_pp"]=100*ratio.loc[x.index]
    x["log_obs"]=np.log1p(x.all_accepted_observation_uuids)
    x["log_species"]=np.log1p(x.all_accepted_recorded_species)
    # Any repeated year of one geographic site is NOT a second distinct
    # within-island site. Require contemporaneous independent geographic
    # lights within at least one true island x year.
    s=x.groupby(["island","year"],observed=True).agg(
         n_sites=("latbin","size"), light_min=(LIGHT,"min"),light_max=(LIGHT,"max"))
    mixed=s.loc[s.n_sites.ge(2) &
                (s.light_max-s.light_min).ge(light_range_min)]
    keep=set(mixed.reset_index().island)
    x=x.loc[x.island.isin(keep)].copy()
    if x.empty:
        # Zero supported same-year within-island light contrasts is
        # an honest NON-ESTIMABLE result, not a failed pandas groupby.
        return x
    n_geo=x.groupby("island",observed=True).apply(
        lambda d:d[["latbin","lonbin"]].drop_duplicates().shape[0],
        include_groups=False)
    x=x.loc[x.island.isin(n_geo.loc[n_geo>=2].index)].copy()
    return x


def weighted_two_way_demean(values,islands,years,weights,tol=1.e-9,max_iter=600):
    """Sparse alternating weighted island/year fixed effects projection."""
    a=np.array(values,dtype=float,copy=True)
    require(a.ndim==2 and len(a)==len(weights),"Wrong residual matrix")
    w=np.asarray(weights,dtype=float)
    require(np.isfinite(w).all() and (w>0).all(),"Invalid source weights")
    ids=[]
    for g in (islands,years):
        codes=pd.factorize(g,sort=True)[0]
        require((codes>=0).all(),"Missing source FE labels")
        ids.append(codes)
    converged=False
    for it in range(max_iter):
        for codes in ids:
            total=np.bincount(codes,weights=w)
            require(np.all(total>0),"Empty fixed effects category")
            mean=np.column_stack([
                np.bincount(codes,weights=w*a[:,j],minlength=len(total))/total
                for j in range(a.shape[1])])
            a-=mean[codes]
        # Numerical criterion is max remaining weighted mean for BOTH FE.
        max_mean=0.
        for codes in ids:
            total=np.bincount(codes,weights=w)
            for j in range(a.shape[1]):
                mean=np.bincount(codes,weights=w*a[:,j],minlength=len(total))/total
                max_mean=max(max_mean,float(np.max(np.abs(mean))))
        if max_mean<tol:
            converged=True
            break
    require(converged,"Two-way weighted FE did not converge")
    return a,{"iterations":it+1,"max_weighted_FE_mean":max_mean}


def fit(frame,weighting="island_equal",adjust_species=False):
    if len(frame)<100 or frame.island.nunique()<30 or frame.year.nunique()<2:
        return {"status":"insufficient_source_support",
                "n_photo_site_years":len(frame),
                "n_islands":int(frame.island.nunique())}
    x=frame.copy()
    if weighting=="island_equal":
        weights=1./x.groupby("island",observed=True).island.transform("size").to_numpy()
    elif weighting=="site_equal":
        weights=np.ones(len(x),dtype=float)
    elif weighting=="accepted_species_weighted":
        weights=x.all_accepted_recorded_species.to_numpy(float)
    else:
        raise ValueError("Unexpected predeclared weighting")
    cols=[LIGHT,"log_obs"]+(["log_species"] if adjust_species else [])
    y=x.frac_pp.to_numpy(float)
    mat=np.column_stack([y,x[cols].to_numpy(float)])
    d, audit=weighted_two_way_demean(mat,x.island.to_numpy(),
                                     x.year.to_numpy(),weights)
    yy=d[:,0]; xx=d[:,1:]; k=len(cols); n=len(x)
    bread_matrix=(xx*weights[:,None]).T@xx
    if np.linalg.matrix_rank(bread_matrix)<k:
        return {"status":"design_rank_deficient","n_islands":x.island.nunique(),
                "n_photo_site_years":n,"predictors":cols}
    bread=np.linalg.inv(bread_matrix)
    beta=bread@((xx*weights[:,None]).T@yy)
    residual=yy-xx@beta
    gids=pd.factorize(x.island,sort=True)[0]
    G=int(gids.max()+1)
    scores=np.zeros((G,k))
    np.add.at(scores,gids,xx*(weights*residual)[:,None])
    cov=bread@(scores.T@scores)@bread
    # Conservative finite-cluster adjustment. Tests remain exploratory.
    n_effective_params=k+G+int(x.year.nunique())-1
    require(n>n_effective_params,"Too few rows after island/year FE")
    cov*=G/(G-1.)*(n-1.)/(n-n_effective_params)
    errs=np.sqrt(np.maximum(np.diag(cov),0))
    q=float(t.ppf(.975,df=G-1))
    results={}
    for j,key in enumerate(cols):
        b=float(beta[j]); se=float(errs[j])
        results[key]={
            "beta_pp":b,"cluster_island_se":se,
            "95CI_pp":[b-q*se,b+q*se],
            "nominal_two_sided_p":float(2*t.sf(abs(b/se),df=G-1)) if se>0 else None}
    return {
      "status":"estimated_exploratory_observational",
      "weights":weighting,"extra_log_accepted_species_adjustment":bool(adjust_species),
      "n_islands":G,"n_site_years":n,"n_source_observation_UUID_sum":int(x.all_accepted_observation_uuids.sum()),
      "year_count":int(x.year.nunique()),"years":sorted(x.year.unique().astype(int).tolist()),
      "island_FE":True,"year_FE":True,"island_cluster_covariance":True,
      "demean_audit":audit,"terms":results,
      "response":"percent of SAME-source WCVP accepted recorded species with GloBI moth-positive floral interaction",
      "NOT_moth_dependence_nor_local_population_decline":True,
      "unmeasured_local_urban_landcover_horticulture_and_GloBI_reporting":True
    }


def run(joined,out):
    out=Path(out)
    require(not out.exists(),"Source-frozen macro result must not overwrite")
    data=pd.read_csv(joined)
    require(REQUIRED.issubset(data.columns),
            "NO true source matched fullflora accepted species denominator")
    primary=prepare(data,obs_min=5,species_min=3)
    variants={
       "island_equal_primary":fit(primary,"island_equal"),
       "site_equal_sensitivity":fit(primary,"site_equal"),
       "island_equal_plus_species_richness":fit(primary,"island_equal",True),
    }
    if len(primary):
        largest=primary.island.value_counts().index[0]
        variants["leave_largest_site_island_out"]=fit(
            primary.loc[~primary.island.eq(largest)],"island_equal")
        variants["strict_10obs_5species"]=fit(
            prepare(data,obs_min=10,species_min=5),"island_equal")
    result={"schema":SCHEMA,"source":"SHA-audited full93m WCVP/GloBI and Zenodo PCNL 1992-2024",
            "years":list(YEARS),
            "source_rows_before_support":len(data),
            "source_islands_before_support":int(data.island.nunique()),
            "retained_photographed_site_years":len(primary),
            "retained_islands_with_sameyear_withinisland_light_contrasts":int(primary.island.nunique()),
            "models":variants,
            "colour_ROI_not_used":True,
            "GloBI_no_record_not_true_nonmoth":True,
            "moth_visit_record_not_effective_or_obligate_moth_pollination":True,
            "historic_PCNL_is_crosscalibrated_past_year_not_local_ground_lux":True,
            "no_plant_local_extinction_or_causal_effect":True,
            "nominal_p_values_post_exploration_not_confirmatory":True}
    out.mkdir(parents=True)
    (out/"ALAN_full93m_moth_GloBI_prior_PCNL_two_way_FE_result.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ALAN_FULL93M_MOTH_POSITIVE_PLANT_SOURCE_PCNL_MACRO_RESULT",
          json.dumps(result,sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joined",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.joined,a.out)

if __name__=="__main__":
    main()
