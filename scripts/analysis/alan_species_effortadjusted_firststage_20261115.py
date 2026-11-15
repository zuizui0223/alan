#!/usr/bin/env python3
"""Exploratory species-specific photographic-record response to historical ALAN.

Scientific scope:
* species recorded vs NOT PHOTOGRAPHICALLY RECORDED among other photographed
  accepted plant taxa, not independently verified biological absence.
* same-island x observation-year photographic sites with real preceding-year
  light contrasts, species-informative strata only, island-equal weights.
* partial out site-level photographic observation effort, with island-year FE.
* population and climate absent from this source: no light-specific causal claim.
* original private species-location and direct botanical trait rows never
  exported; only aggregated diagnostics and group effects are written.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t as student_t, spearmanr

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, load_sites, load_labels, original_accepted_species, require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import (
    direct_plant_traits
)

SPANS = (.02, .10)
EXPECTED = dict(source_islands=562, source_cells=77647,
                source_events=2012925, source_species=14330,
                common_islands=243, common_cells=72472,
                strict_selfing_species=703)
MIN_CLUSTER_SE=10
MIN_CLUSTER_ASSESS=30

def island_fold(island):
    return hashlib.sha256(str(island).encode()).digest()[0]&1

def prepare_support(sites, flora, traits, *, span=.02):
    """True island-year photo-gradient screen with no biology-as-absence."""
    require(not sites.duplicated(KEY).any(),"Duplicated source photo-cell")
    require(not flora.duplicated(KEY+["species"]).any(),"Duplicated source species-site")
    src=sites.copy().reset_index(drop=True)
    src["sid"]=np.arange(len(src),dtype=np.int32)
    g=(src.groupby(["island","year"],observed=True)[LIGHT]
       .agg(n_sites="size",minimum="min",maximum="max",median="median")
       .reset_index())
    g=g.loc[g.n_sites.ge(6)&(g.maximum-g.minimum).ge(span)]
    site=src.merge(g[["island","year","median"]],
                   on=["island","year"],validate="many_to_one")
    site["hi"]=(site[LIGHT]>site["median"]).astype("int8")
    halves=site.groupby(["island","year","hi"],observed=True).size().unstack(
        "hi",fill_value=0)
    for k in (0,1):
        if k not in halves:halves[k]=0
    good=halves.loc[halves[0].ge(2)&halves[1].ge(2)].reset_index()[
        ["island","year"]]
    site=site.merge(good,on=["island","year"],validate="many_to_one")
    site["iy"],_=pd.factorize(site.island.astype(str).str.cat(
        site.year.astype(str),sep="|"),sort=True)
    site["fold"]=site.island.map(island_fold).astype("int8")
    rows=flora[KEY+["species"]].merge(
        site[KEY+["sid","iy","hi","fold"]],
        on=KEY,how="inner",validate="many_to_one")
    rows["lo"]=1-rows["hi"]
    nstr=site.groupby("iy",observed=True).size().rename("stratum_sites")
    counts=(rows.groupby(["species","iy","island","year","fold"],
                         observed=True)
            .agg(n=("species","size"),lo=("lo","sum"),hi=("hi","sum"))
            .reset_index())
    counts["stratum_sites"]=counts.iy.map(nstr)
    counts=counts.loc[counts.n.lt(counts.stratum_sites)].copy()
    require((counts.n.gt(0)&counts.stratum_sites.gt(counts.n)).all(),
            "Noninformative species-strata admitted")
    summary=(counts.groupby("species",observed=True)
             .agg(n=("n","sum"),n_strata=("n","size"),
                  n_islands=("island","nunique"),n_years=("year","nunique"),
                  positive_low=("lo","sum"),positive_high=("hi","sum")))
    folds=(counts.groupby(["species","fold"],observed=True)
           .agg(n=("n","sum"),n_strata=("n","size")).reset_index())
    folds["ok"]=folds.n.ge(10)&folds.n_strata.ge(2)
    summary["both_fold"]=folds.groupby("species").ok.sum().reindex(
        summary.index,fill_value=0).ge(2)
    summary["strict"]=(summary.n.ge(30)&summary.n_strata.ge(6)&
                      summary.n_islands.ge(2)&summary.n_years.ge(2)&
                      summary.positive_low.ge(5)&summary.positive_high.ge(5)&
                      summary.both_fold)
    selected=summary.loc[summary.strict].index.intersection(
        traits.loc[traits.autonomy.isin(["autonomous","absent"]),"species"])
    selected=set(selected)
    pos=rows.loc[rows.species.isin(selected),["species","sid","iy","fold"]]
    pos=pos.sort_values(["species","sid"])
    eligible_counts=counts.loc[counts.species.isin(selected)]
    eligible_counts=eligible_counts.sort_values(["species","iy"])
    return site,pos,eligible_counts,selected,{
        "gradient_span_min_log1p":span,
        "eligible_islands":int(site.island.nunique()),
        "eligible_photo_site_years":len(site),
        "island_years":int(len(good)),
        "source_photographic_plant_species":int(rows.species.nunique()),
        "strict_selfing_species":len(selected),
        "strict_selfing_autonomous":int(
            traits.loc[traits.species.isin(selected),"autonomy"].eq("autonomous").sum()),
        "strict_selfing_absent":int(
            traits.loc[traits.species.isin(selected),"autonomy"].eq("absent").sum())
    }

def site_arrays(site):
    z=site.copy().reset_index(drop=True)
    # Source has no independent population density / climate covariates.
    for c in ("population_density","worldpop","temperature","precipitation"):
        require(c not in z.columns,"Source covariate schema changed; model must be refrozen")
    z["effort"]=z.log_obs.astype(float)
    z["x"]=z[LIGHT].astype(float)
    means=z.groupby("iy",observed=True)[["x","effort"]].transform("mean")
    z["xc"]=z.x-means.x
    z["ec"]=z.effort-means.effort
    require(np.isfinite(z[["x","effort","xc","ec"]]).all().all(),
            "Nonfinite exposures or effort")
    # True same-stratum partial light~effort correlation; distinct island-year FE.
    corr=float(np.corrcoef(z.xc,z.ec)[0,1])
    return z,{"pcnl_effort_within_island_year_correlation":corr,
              "pcnl_effort_partial_r2":corr**2,
              "population_density_covariate":"NOT_PRESENT_IN_AUTHENTIC_SOURCE",
              "climate_covariates":"NOT_PRESENT_IN_AUTHENTIC_SOURCE",
              "response":"photographically_observed_species_vs_photographic_nonrecord",
              "light":"rolling_mean_previous_calendar_year_log1p_PCNL",
              "effort":"log1p_accepted_photographic_observation_UUIDs",
              "island_year_fixed_effects":True}

def fit_one(site, full_strata, pos_site_ids, *, fold=None):
    """2-regressor weighted within-island-year LPM, island-cluster sandwich.

    Include full photographed survey-opportunity sites for each informative
    species stratum; exclude whole strata with species never photographed.
    """
    if fold is not None:
        s=site.loc[site.fold.eq(fold)&site.iy.isin(full_strata)].copy()
    else:
        s=site.loc[site.iy.isin(full_strata)].copy()
    if len(s)<20:return {"ok":False,"why":"too_few_site_opportunities"}
    idx=s.sid.to_numpy(np.int64)
    y=np.isin(idx,pos_site_ids,assume_unique=True).astype(float)
    means=s.assign(y=y).groupby("iy",observed=True).y.transform("mean").to_numpy()
    yc=y-means
    x=s.xc.to_numpy(float)
    e=s.ec.to_numpy(float)
    islands,codes=pd.factorize(s.island,sort=True)
    # island-equal TOTAL sampling weight within the species-specific sample
    # pd.factorize returns (codes, uniques), not vice versa
    codes,levels=pd.factorize(s.island,sort=True)
    G=len(levels)
    w=1/np.bincount(codes).astype(float)[codes]
    X=np.column_stack([x,e])
    xtwx=(X.T*w)@X
    if not np.isfinite(xtwx).all() or np.linalg.matrix_rank(xtwx)<2:
        return {"ok":False,"why":"rank_deficient_light_effort"}
    condition=float(np.linalg.cond(xtwx))
    if condition>1.e10:
        return {"ok":False,"why":"ill_conditioned_light_effort"}
    bread=np.linalg.inv(xtwx)
    beta=bread@(X.T@(w*yc))
    resid=yc-X@beta
    if not np.isfinite(beta).all():
        return {"ok":False,"why":"nonfinite_beta"}
    result={"ok":True,"beta_pp":float(100*beta[0]),
            "effort_beta_pp":float(100*beta[1]),
            "n_islands":G,"n_strata":int(s.iy.nunique()),
            "n_source_photo_site_years":int(len(s)),
            "n_species_recorded_site_years":int(y.sum()),
            "species_recorded_fraction":float(y.mean()),
            "condition":condition}
    if G<MIN_CLUSTER_SE:
        result.update(cluster_se_status="INSUFFICIENT_ISLAND_CLUSTER",
                      se_pp=None,ci95_pp=None,p_nominal=None)
        return result
    score=np.zeros((G,2),dtype=float)
    np.add.at(score,codes,X*(w*resid)[:,None])
    degfree=len(s)-s.iy.nunique()-2
    if degfree<=0:
        result.update(cluster_se_status="NONPOSITIVE_DF",
                      se_pp=None,ci95_pp=None,p_nominal=None)
        return result
    cov=((G/(G-1))*((len(s)-1)/degfree))*bread@(score.T@score)@bread
    se=float(np.sqrt(max(cov[0,0],0)))
    tcrit=float(student_t.ppf(.975,G-1))
    result.update(cluster_se_status="ISLAND_CR1_EXPLORATORY",
                  se_pp=100*se,
                  ci95_pp=[float(100*(beta[0]-tcrit*se)),
                           float(100*(beta[0]+tcrit*se))],
                  p_nominal=float(2*student_t.sf(abs(beta[0]/se),G-1))
                  if se>0 else None)
    return result

def percentile(v):
    a=np.asarray(list(v),float)
    a=a[np.isfinite(a)]
    if not len(a):return {"n":0}
    return {"n":len(a),"mean":float(np.mean(a)),
            "median":float(np.median(a)),
            "q025":float(np.quantile(a,.025)),
            "q25":float(np.quantile(a,.25)),
            "q75":float(np.quantile(a,.75)),
            "q975":float(np.quantile(a,.975))}

def group_summary(fits,group):
    ok=fits.loc[fits.ok & fits.autonomy.eq(group)].copy()
    reliable=ok.loc[ok.n_islands.ge(MIN_CLUSTER_ASSESS)&ok.se_pp.notna()].copy()
    return {
        "n_total_tagged":int(fits.autonomy.eq(group).sum()),
        "n_slope_fitted":len(ok),
        "n_with_ge10_island_cluster_se":int((ok.se_pp.notna()).sum()),
        "n_with_ge30_islands":len(reliable),
        "unweighted_slope_pp":percentile(ok.beta_pp),
        "ge30_islands_slope_pp":percentile(reliable.beta_pp),
        "n_nominal_ci_below_0_ge30_islands":int(reliable.ci95_pp.map(
            lambda ci:ci is not None and ci[1]<0).sum()),
        "n_nominal_ci_above_0_ge30_islands":int(reliable.ci95_pp.map(
            lambda ci:ci is not None and ci[0]>0).sum()),
    }

def genus_adjusted_difference(fits,min_islands=0):
    d=fits.loc[fits.ok].copy()
    if min_islands:d=d.loc[d.n_islands.ge(min_islands)]
    if d.empty:return {"status":"NOT_ESTIMABLE_NO_SPECIES"}
    mixed=(d.groupby("genus",observed=True).autonomy.nunique())
    valid=set(mixed[mixed.eq(2)].index)
    d=d.loc[d.genus.isin(valid)].copy()
    if not len(d):return {"status":"NOT_ESTIMABLE_NO_GENUS_WITH_BOTH_SELFING_STATES"}
    y=d.beta_pp.to_numpy(float)
    x=d.autonomy.eq("autonomous").astype(float).to_numpy()
    means=d.assign(x=x).groupby("genus")[["x","beta_pp"]].transform("mean")
    xx=x-means.x.to_numpy()
    yy=y-means.beta_pp.to_numpy()
    den=float(xx@xx)
    if den<1.e-8:return {"status":"NOT_ESTIMABLE_WITHIN_GENUS_NO_VARIATION"}
    return {"status":"EXPLORATORY_DESCRIPTIVE_GENUS_FIXED_EFFECT",
            "beta_autonomous_minus_absent_pp":float(xx@yy/den),
            "n_mixed_genera":len(valid),"n_species":len(d),
            "no_valid_causal_or_clustered_inference":True}

def run(sites_path,globi_path,groups,ledger,out,span=.02):
    require(not out.exists(),"Refuse overwrite")
    sites=load_sites(sites_path)
    flora=original_accepted_species(groups,sites,load_labels(globi_path))
    require((sites.island.nunique(),len(sites),len(flora),
             flora.species.nunique())==(562,77647,2012925,14330),
            "Original source census changed")
    traits=direct_plant_traits(ledger,high_only=False)
    high=direct_plant_traits(ledger,high_only=True)
    taxa=pd.DataFrame({"species":flora.species.drop_duplicates()})
    taxa["autonomy"]=taxa.species.map(traits.autonomous_selfing_capacity)
    taxa["high_quality"]=taxa.species.map(high.autonomous_selfing_capacity).isin(
        ["autonomous","absent"])
    taxa["genus"]=taxa.species.str.split().str[0]
    site,positive,groupspec,selected,support=prepare_support(
        sites,flora,taxa,span=span)
    if span==.02:
        require((support["eligible_islands"],
                support["eligible_photo_site_years"],
                support["strict_selfing_species"],
                support["strict_selfing_autonomous"],
                support["strict_selfing_absent"])==
                (243,72472,703,489,214),"Frozen source selection changed")
    site,geo=site_arrays(site)
    # Group extraction is restricted to pre-screened species.
    pos={sp:np.asarray(g.sid,dtype=np.int64) for sp,g in
         positive.groupby("species",observed=True,sort=True)}
    info={sp:np.asarray(g.iy,dtype=int) for sp,g in
          groupspec.groupby("species",observed=True,sort=True)}
    traits_subset=taxa.loc[taxa.species.isin(selected)].set_index("species")
    fits=[]
    for j,sp in enumerate(sorted(selected)):
        st=info[sp]
        f=fit_one(site,st,pos[sp])
        fold0=fit_one(site,st,pos[sp],fold=0)
        fold1=fit_one(site,st,pos[sp],fold=1)
        if f["ok"]:
            record={"species":sp,"genus":str(traits_subset.loc[sp,"genus"]),
                    "autonomy":str(traits_subset.loc[sp,"autonomy"]),
                    "high_quality":bool(traits_subset.loc[sp,"high_quality"]),
                    "ok":True,**f,
                    "fold0_beta_pp":fold0.get("beta_pp"),
                    "fold1_beta_pp":fold1.get("beta_pp"),
                    "fold0_islands":fold0.get("n_islands",0),
                    "fold1_islands":fold1.get("n_islands",0)}
        else:
            record={"species":sp,"genus":str(traits_subset.loc[sp,"genus"]),
                    "autonomy":str(traits_subset.loc[sp,"autonomy"]),
                    "high_quality":bool(traits_subset.loc[sp,"high_quality"]),
                    "ok":False,"why":f.get("why")}
        fits.append(record)
        if (j+1)%100==0:
            print("SOURCE_SPECIES_FIT_PROGRESS",j+1,"of",len(selected),flush=True)
    df=pd.DataFrame(fits)
    ok=df.loc[df.ok].copy()
    if not ok.empty:
        for f in ("fold0_beta_pp","fold1_beta_pp"):
            ok[f]=pd.to_numeric(ok[f],errors="coerce")
        folds=ok.dropna(subset=["fold0_beta_pp","fold1_beta_pp"])
        r=float(spearmanr(folds.fold0_beta_pp,folds.fold1_beta_pp).statistic) if len(folds)>2 else None
        sign=((np.sign(folds.fold0_beta_pp)==np.sign(folds.fold1_beta_pp)).mean()
              if len(folds)>0 else None)
        fold_a=folds.loc[folds.autonomy.eq("autonomous")]
        fold_b=folds.loc[folds.autonomy.eq("absent")]
        meta={
            "n_valid_fold_beta":len(folds),
            "fold_slope_sign_agreement_fraction":float(sign) if sign is not None else None,
            "fold_slope_rank_spearman":r,
            "autonomous_full_unweighted_mean_pp":float(
                ok.loc[ok.autonomy.eq("autonomous"),"beta_pp"].mean()),
            "absent_full_unweighted_mean_pp":float(
                ok.loc[ok.autonomy.eq("absent"),"beta_pp"].mean()),
            "autonomous_minus_absent_unweighted_mean_pp":float(
                ok.loc[ok.autonomy.eq("autonomous"),"beta_pp"].mean()-
                ok.loc[ok.autonomy.eq("absent"),"beta_pp"].mean()),
            "autonomous_minus_absent_fold0_mean_pp":(
                float(fold_a.fold0_beta_pp.mean()-fold_b.fold0_beta_pp.mean())
                if len(fold_a) and len(fold_b) else None),
            "autonomous_minus_absent_fold1_mean_pp":(
                float(fold_a.fold1_beta_pp.mean()-fold_b.fold1_beta_pp.mean())
                if len(fold_a) and len(fold_b) else None),
            "genus_FE_all":genus_adjusted_difference(ok),
            "genus_FE_30plus_islands":genus_adjusted_difference(ok,min_islands=30),
            "independent_beta_meta_se_NOT_computed":True,
            "country_geography_block_bootstrap_NOT_computed":True,
        }
        onlyhigh=ok.loc[ok.high_quality]
        meta["high_quality_autonomous_minus_absent_unweighted_pp"]=(
            float(onlyhigh.loc[onlyhigh.autonomy.eq("autonomous"),"beta_pp"].mean()-
                  onlyhigh.loc[onlyhigh.autonomy.eq("absent"),"beta_pp"].mean())
            if len(onlyhigh) else None)
    else:
        meta={"status":"NO_SLOPES_ESTIMABLE"}
    summary={
        "schema":"ALAN_genuine_island_species_LPM_firststage_effort_and_islandyear_v1",
        "status":"ORIGINAL_SOURCE_EXPLORATORY_SPECIES_FITS_NO_CAUSAL_OR_FULL_META_CLAIM",
        "base_original_source":{"islands":562,"photo_site_years":77647,
           "species_site_years":2012925,"species":14330,
           "direct_species_trait_global_sha256":
               "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"},
        "support":support,"site_covariate_audit":geo,
        "fitting":{
          "n_selected":len(df),"n_successful":len(ok),
          "n_nonfit":int((~df.ok).sum()),
          "n_with_ge10_island_cluster_SE":int(ok.se_pp.notna().sum()) if len(ok) else 0,
          "n_with_ge30_islands":int(ok.n_islands.ge(30).sum()) if len(ok) else 0,
          "unfittable_reason_counts":df.loc[~df.ok,"why"].value_counts().to_dict()
               if (~df.ok).any() else {},
          "all":group_summary(df,"autonomous"),
          "absent":group_summary(df,"absent"),
          "high_quality_successes":int(ok.high_quality.sum()) if len(ok) else 0,
          "all_species_slope_pp":percentile(ok.beta_pp) if len(ok) else {"n":0}
        },
        "descriptive_second_stage":meta,
        "limits":{
          "no_true_species_absence":True,
          "population_density_and_climate_unavailable":True,
          "urbanization_and_plant_habitat_confounding_unresolved":True,
          "no_ground_lux_or_manipulated_ALAN":True,
          "species_slope_p_values_multiple_unadjusted":True,
          "not_global_mainland_comparison":True,
          "not_pollination_or_seed_reproduction_response":True,
          "no_independent_phylogenetic_or_shared_site_bootstrap_uncertainty":True,
          "selection_on_photo_observation_and_documented_trait":True,
          "partial_fit_stage_only_not_claiming_meta_confirmed":True
        }
    }
    out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_FIRSTSTAGE_AGGREGATES",json.dumps({
        "support":support,"fitting":summary["fitting"],
        "descriptive_second_stage":meta,"covariate":geo
    },ensure_ascii=False,allow_nan=False),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("sites","globi","groups","ledger","out"):
        p.add_argument("--"+a,required=True,type=Path)
    p.add_argument("--span",type=float,default=.02)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.out,a.span)
