#!/usr/bin/env python3
"""Worldwide 2013-2025 island GBIF 4-family SOURCE COMPOSITION vs PRIOR PCNL.

Primary outcome: within-original-site/year PHOTO-RECORD fraction of explicitly
identified moth families Noctuidae+Geometridae among those two and Nymphalidae+
Pieridae; this is NOT real moth abundance, day/night activity, pollination,
wild plant persistence or unbiased community richness.

Observational diagnostics only. Exact original 2026-08 GSHHG island source
snapshot; island×observation-year fixed intercept and island-cluster robust
nominal uncertainty; species absence never imputed from GBIF source negatives.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from scripts.analysis.alan_global_GBIF_true4family_priorPCNL_20261110 import (
    SCHEMA as LIGHT_SCHEMA,sha256,require)
SOURCE_ROWS=228361
SOURCE_RECORDS=22713063
MIN_OBSERVATION_OPPORTUNITY=5
MIN_SAME_YEAR_LIGHT_SPAN=.02
SCHEMA="ALAN_worldwide_island_GBIF_fourfamily_prior_PCNL_photo_record_fraction_FE_v1"

def source(source,receipt):
    info=json.loads(Path(receipt).read_text())
    require(info.get("schema")==LIGHT_SCHEMA and
            info["true_global_original_source_cells"]==SOURCE_ROWS and
            info["true_global_fourfamily_original_human_observation_records"]==SOURCE_RECORDS and
            sha256(source)==info["real_complete_joined_global_source_CSV_SHA256"],
            "Not genuine prior-year-PCNL joined 2026-08 full GSHHG source")
    x=pd.read_csv(source,dtype={"island":str})
    k=("island","year","latbin","lonbin")
    require(len(x)==SOURCE_ROWS and not x.duplicated(list(k)).any() and
            (x.PCNL_year==x.year-1).all() and
            int(x.all_four_family_records.sum())==SOURCE_RECORDS,
            "Globally photographed observation/year source altered")
    require((x.moth_family_records+x.butterfly_family_records==
             x.all_four_family_records).all() and
            (x.all_four_family_records>0).all(),
            "Original GBIF documented Lepidoptera categories non-additive")
    x["recorded_moth_family_share_pp"]=100*x.moth_family_records/x.all_four_family_records
    x["log_all_fourfamily_observations"]=np.log1p(x.all_four_family_records)
    x["island_year"]=x.island+"|"+x.year.astype(str)
    x["centre_lat"]=.05*(x.latbin.to_numpy(float)+.5)-90.
    x["log10_source_island_area"]=np.log10(x.area_km2)
    require(x.recorded_moth_family_share_pp.between(0,100).all() and
            x.log10_source_island_area.notna().all(),
            "Recorded composition outside valid source range")
    return x,info

def frame(d,response,include_zero_category=False,obs_threshold=5):
    z=d.loc[d.prior_log1p_PCNL.notna() &
            d.all_four_family_records.ge(obs_threshold)].copy()
    if not include_zero_category:
        z=z.loc[z.moth_family_records.gt(0)&z.butterfly_family_records.gt(0)].copy()
    z=z.loc[z[response].notna() &
            z.centre_lat.notna() &
            z.area_km2.gt(0) &
            z.prior_log1p_PCNL.notna()].copy()
    support=z.groupby("island_year",observed=True).agg(
        true_distinct_source_cells=("lonbin","size"),
        n_light_levels=("prior_log1p_PCNL","nunique"),
        first_light=("prior_log1p_PCNL","min"),
        last_light=("prior_log1p_PCNL","max"))
    good=support.index[(support.true_distinct_source_cells>=2)&
         (support.n_light_levels>=2)&
         ((support.last_light-support.first_light)>=MIN_SAME_YEAR_LIGHT_SPAN)]
    z=z.loc[z.island_year.isin(good)].copy()
    info={
        "supported_independent_islands":int(z.island.nunique()),
        "supported_original_site_years":int(len(z)),
        "same_island_actual_year_groups_with_2plus_geo_sites":int(len(good)),
        "source_cell_year_before_light_positivity":int(
             d.prior_log1p_PCNL.notna().sum()),
        "source_min_guild_observation_records":obs_threshold,
        "source_conditional_on_both_recorded_families":not include_zero_category,
        "year_range":sorted(z.year.unique().astype(int).tolist()),
        "no_missing_moth_source_claim_treated_as_biological_absence":True,
    }
    return z,info

def transform_weights(x,mode):
    if mode=="island_equal":
        w=1./x.groupby("island").island.transform("size").to_numpy(float)
    elif mode=="site_equal":
        w=np.ones(len(x),dtype=float)
    elif mode=="island_year_equal":
        w=1./x.groupby("island_year").island_year.transform("size").to_numpy(float)
    else:
        raise ValueError("Unregistered photo-source cell weighting")
    require(np.isfinite(w).all() and (w>0).all(),"Weight not source defined")
    return w

def original_group_demean(values,labels,w):
    code=pd.factorize(labels,sort=True)[0]
    require(np.all(code>=0) and values.shape[0]==len(w),
            "Unknown island×year fixed effect key")
    totals=np.bincount(code,weights=w)
    require(np.all(totals>0),"Island×year group has zero sampled cells")
    means=np.column_stack([
        np.bincount(code,weights=w*values[:,j],minlength=len(totals))/totals
        for j in range(values.shape[1])])
    out=values-means[code]
    check=np.column_stack([
       np.bincount(code,weights=w*out[:,j],minlength=len(totals))/totals
       for j in range(out.shape[1])])
    require(np.max(np.abs(check))<1e-8,
            "Weighted original island×year fixed effects not removed")
    return out

def fit(d,response,weighting="island_equal",adjust_obs=True,
        include_zero_category=False,include_area_modifier=False,
        obs_threshold=5,remove_largest_island=None):
    x,support=frame(d,response,include_zero_category,obs_threshold)
    if remove_largest_island is not None:
        x=x.loc[~x.island.eq(remove_largest_island)].copy()
        support["removed_one_source_popular_island"]=str(remove_largest_island)
    G=x.island.nunique();N=len(x);H=x.island_year.nunique()
    if G<30 or N<max(180,H+15):
        return {"status":"insufficient_independent_real_islands_or_source_cells",
                "support":support,"status_counts":{"independent_islands":int(G),
                    "eligible_site_years":N,"true_island_year_groups":int(H)}}
    x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
    cols=["prior_log1p_PCNL","centre_lat_c"]
    if adjust_obs:cols.append("log_all_fourfamily_observations")
    if include_area_modifier:
        unit_area=x.groupby("island").log10_source_island_area.first()
        center=float(unit_area.median())
        x["light_by_log10_area_c"]=x.prior_log1p_PCNL*(x.log10_source_island_area-center)
        cols.append("light_by_log10_area_c")
    y=x[response].to_numpy(float)
    X=x[cols].to_numpy(float)
    w=transform_weights(x,weighting)
    data=original_group_demean(np.column_stack([y,X]),x.island_year.to_numpy(),w)
    yy=data[:,0];xx=data[:,1:]
    k=len(cols)
    bread_inv=(xx*w[:,None]).T@xx
    if np.linalg.matrix_rank(bread_inv)<k:
        return {"status":"source_design_rank_deficient",
                "predictors":cols,"support":support,
                "n_independent_islands":int(G)}
    bread=np.linalg.inv(bread_inv)
    b=bread@((xx*w[:,None]).T@yy)
    residual=yy-xx@b
    labels=pd.factorize(x.island,sort=True)[0]
    score=np.zeros((G,k))
    np.add.at(score,labels,xx*(w*residual)[:,None])
    cov=bread@(score.T@score)@bread
    require(N>H+k and G>1,
            "Source cluster covariance degrees of freedom unidentified")
    cov*=G/(G-1)*(N-1)/(N-H-k)
    se=np.sqrt(np.maximum(0,np.diag(cov)))
    q=float(t.ppf(.975,G-1))
    params={}
    for i,col in enumerate(cols):
        params[col]={
            "beta":float(b[i]),"cluster_island_SE":float(se[i]),
            "nominal_95CI":[float(b[i]-q*se[i]),float(b[i]+q*se[i])],
            "nominal_two_sided_unadjusted_p":float(2*t.sf(abs(b[i]/se[i]),G-1))
                if se[i]>0 else None}
    R2=1-(sum(w*residual**2)/sum(w*yy**2))
    support.update({"observed_true_source_islands_in_fit":int(G),
        "observed_true_island_year_fixed_effects":int(H),
        "true_siteyear_rows_after_source_controls":N})
    return {
        "status":"estimated_exploratory_observation_composition_only",
        "response":response,
        "main_exposure":"log1p(annual PCNL source satellite pixel year-before recorded insects)",
        "true_island_year_fixed_intercepts":True,
        "actual_GSHHG_island_cluster_covariance":True,
        "weights":weighting,
        "adjusting_actual_recorded_four_family_source_opportunity":adjust_obs,
        "interaction_with_source_island_log10_area":include_area_modifier,
        "sampling_supported_by_both_moth_and_butterfly_record_categories":
           not include_zero_category,
        "original_source_sample_size":N,"observed_islands":int(G),
        "island_year_FEs":int(H),
        "source_support":support,
        "within_island_year_weighted_R2":float(R2),
        "predictor_coefficients":params,
        "NOT_causal_moth_abundance_or_effective_night_pollination":True,
        "independent_recording_effort_and_horticulture_unmeasured":True}

def run(source_csv,receipt,out):
    all_source,source_provenance=source(source_csv,receipt)
    primary="recorded_moth_family_share_pp"
    main,info=frame(all_source,primary,False,5)
    if len(main):
        largest_island=main.island.value_counts().index[0]
    else:
        largest_island=None
    specs={
        "island_equal_primary_recorded_guild_fraction":
            fit(all_source,primary,"island_equal",True,False,False,5),
        "site_equal_weighting":
            fit(all_source,primary,"site_equal",True,False,False,5),
        "island_year_equal_weighting":
            fit(all_source,primary,"island_year_equal",True,False,False,5),
        "extra_island_area_moderation":
            fit(all_source,primary,"island_equal",True,False,True,5),
        "source_effort_not_adjusted":
            fit(all_source,primary,"island_equal",False,False,False,5),
        "include_recorded_zero_category_but_not_true_absence":
            fit(all_source,primary,"island_equal",True,True,False,5),
        "strict_20_source_records_and_both":
            fit(all_source,primary,"island_equal",True,False,False,20),
        "all_family_record_count_negative_control":
            fit(all_source,"log_all_fourfamily_observations","island_equal",
                False,False,False,5),
    }
    if largest_island is not None:
        specs["omit_largest_documented_photo_island"]=fit(
            all_source,primary,"island_equal",True,False,False,5,
            remove_largest_island=largest_island)
    report={
        "schema":SCHEMA,
        "source_original_years":list(range(2013,2026)),
        "real_global_original_source_record_count":SOURCE_RECORDS,
        "real_global_original_site_year_cells":SOURCE_ROWS,
        "original_source_GBIF_DOI":"10.15468/dl.8yrbe7",
        "original_historical_PCNL_DOI":"10.5281/zenodo.17013414",
        "source_signed_original_source_CSV_SHA256":sha256(source_csv),
        "original_PCNl_join_source_summary":source_provenance,
        "all_estimates_are_documented_photo_record_composition_not_insect_abundance":True,
        "original_GBIF_moth_families":["Noctuidae","Geometridae"],
        "original_GBIF_butterfly_families":["Nymphalidae","Pieridae"],
        "source_spatial_cell_size":".05 decimal degrees, not field survey 500m",
        "pre_year_exposure_not_contemporaneous_visit_light":True,
        "no_true_absent_moth_or_butterfly_species_inferred":True,
        "same_GBIF_observer_platform_and_detection_controls_unmeasured":True,
        "island_geography_sampling_weights_not_probability_design":True,
        "source_descriptive_postexploratory_nominal_CI_not_confirmatory":True,
        "all_prespecified_source_model_sensitivities":specs}
    out=Path(out)
    require(not out.exists(),"Original source model is append-only")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_GLOBAL_GBIF_MOTH_PHOTO_RECORD_VS_PRIOR_LIGHT",
          json.dumps({"primary":specs["island_equal_primary_recorded_guild_fraction"],
                      "negative_control":specs["all_family_record_count_negative_control"],
                      "area_modifier":specs["extra_island_area_moderation"]},
                     sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();run(a.source,a.receipt,a.out)
