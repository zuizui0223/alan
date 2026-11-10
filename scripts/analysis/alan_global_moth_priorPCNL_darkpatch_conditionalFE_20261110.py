#!/usr/bin/env python3
"""Observational global moth-record source fraction versus historic dark topology.

Actual 2026-08 GBIF moth-v-butterfly PHOTO record composition (not trapping),
actual pre-observation annual 2012-24 PCNL 5km island-land mask, 196 islands
max. Within island×same photo year, adjusting local source PCNL pixel,
source dark fraction, true source 4-connected dark patch fraction,
latitude and recorded source photo opportunity.
No effective nocturnal moth service or population density is inferred.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
from scripts.analysis.alan_global_trueGBIF_PCNl_13year_darkpatch_20261110 import (
    SCHEMA as PATCH_SCHEMA,digest as sha,require,source_eligible)
from scripts.analysis.alan_global_GBIF_moth_record_fraction_within_island_FE_20261110 import (
    original_group_demean,transform_weights)

SCHEMA="ALAN_global_GBIF_recorded_moth_fraction_vs_prior_PCNL_darkpatch_FE_v1"
THRESHOLDS=(0.,3.)
REGISTRY=("prior_log1p_PCNL","dark_land_fraction",
          "largest_dark_component_fraction_of_dark",
          "log_all_fourfamily_observations","centre_lat_c")

def eligible_full(source_csv,source_receipt,patch_csv,patch_receipt):
    rec=json.loads(Path(patch_receipt).read_text())
    require(rec["schema"]==PATCH_SCHEMA and
            rec["source_actual_global_GBIF_photo_cells"]==62688 and
            rec["source_original_islands"]==196 and
            sha(patch_csv)==rec["source_2threshold_global_patch_CSV_SHA256"],
            "No original verified 13year/5km darkpatch source receipt")
    photo=source_eligible(source_csv,source_receipt)
    a=pd.read_csv(patch_csv,dtype={"island":str})
    keys=["island","year","latbin","lonbin"]
    require(len(a)==2*len(photo) and
            set(a.threshold_PCNL_native)==set(THRESHOLDS) and
            not a.duplicated(keys+["threshold_PCNL_native"]).any(),
            "Observed-source dark-patch repeated/photo IDs changed")
    for t in THRESHOLDS:
        s=a.loc[a.threshold_PCNL_native.eq(t)]
        require(len(s)==len(photo) and
                not s.duplicated(keys).any(),
                "Original source threshold lacks genuine photo cell")
    joined=a.merge(photo[keys+["recorded_moth_family_share_pp",
                              "log_all_fourfamily_observations",
                              "prior_log1p_PCNL",
                              "centre_lat","area_km2","island_year"]],
                   on=keys,validate="many_to_one",how="left",indicator=True)
    require(joined._merge.eq("both").all(),
            "Real source biology/photo observation not joined to exact past light patch")
    return joined,rec

def source_stable_group_selection(x,threshold,min_connect_span=.05,require_dark_area_span=False):
    d=x.loc[x.threshold_PCNL_native.eq(threshold) &
            x.status.eq("observed_satellite_geometry_only") &
            x.n_dark_raster_pixels.gt(0) &
            x.largest_dark_component_fraction_of_dark.notna() &
            x.dark_land_fraction.notna() &
            x.prior_log1p_PCNL.notna()].copy()
    d["centre_lat_c"]=d.centre_lat-d.centre_lat.median()
    g=d.groupby("island_year",observed=True).agg(
        n_sites=("lonbin","size"),
        n_light_levels=("prior_log1p_PCNL","nunique"),
        light_span=("prior_log1p_PCNL",lambda z:float(z.max()-z.min())),
        dark_span=("dark_land_fraction",lambda z:float(z.max()-z.min())),
        connection_span=("largest_dark_component_fraction_of_dark",
                         lambda z:float(z.max()-z.min())))
    supported=(g.n_sites>=2)&(g.n_light_levels>=2)&(g.light_span>=.02)&(
        g.connection_span>=min_connect_span)
    if require_dark_area_span:
        supported &= g.dark_span>=.05
    x=d.loc[d.island_year.isin(g.index[supported])].copy()
    audit={"actual_threshold":threshold,
           "full_genuine_darkpatch_rows":len(d),
           "source_islands_with_real_darkpatches":int(d.island.nunique()),
           "retained_islandyear_groups_with_independent_dark_connectivity_contrast":
              int(supported.sum()),
           "retained_actual_independent_islands":int(x.island.nunique()),
           "retained_actual_photo_grid_years":len(x),
           "prespecified_min_islandyear_connected_dark_patch_span":min_connect_span,
           "additional_min_dark_fraction_span_if_requested":require_dark_area_span,
           "no_dark_pixel_sites_are_censored_not_zero_connectivity":True}
    return x,audit

def fit_source(x,y,weighting,columns):
    n=len(x);G=x.island.nunique();H=x.island_year.nunique();k=len(columns)
    if n<max(200,H+10) or G<30:
        return {"status":"insufficient_original_source_islands_for_joint_dark_effect",
                "source_islands":int(G),"source_siteyear_cells":n,
                "source_year_groups":int(H)}
    xv=x[columns].to_numpy(float)
    yv=x[y].to_numpy(float)
    require(np.isfinite(xv).all() and np.isfinite(yv).all(),
            "Original dark-buffers have unresolved missing source values")
    weights=transform_weights(x,weighting)
    z=original_group_demean(np.column_stack([yv,xv]),x.island_year.to_numpy(),weights)
    yy=z[:,0];xx=z[:,1:]
    matrix=(xx*weights[:,None]).T@xx
    cond=float(np.linalg.cond(matrix))
    if np.linalg.matrix_rank(matrix)<k or cond>1e11:
        return {"status":"joint_dark_fraction_and_connectivity_rank_or_condition_failure",
                "n_true_islands":int(G),"n_source_siteyear_cells":n,
                "source_predictor_condition_number":cond}
    bread=np.linalg.inv(matrix)
    b=bread@((xx*weights[:,None]).T@yy)
    err=yy-xx@b
    islands=pd.factorize(x.island,sort=True)[0]
    scores=np.zeros((G,k))
    np.add.at(scores,islands,xx*(weights*err)[:,None])
    cov=bread@(scores.T@scores)@bread
    require(n>H+k and G>1,
            "Missing independent source island cluster covariance degrees of freedom")
    cov*=G/(G-1)*(n-1)/(n-H-k)
    se=np.sqrt(np.maximum(0,np.diag(cov)))
    cutoff=float(t.ppf(.975,G-1))
    terms={}
    for i,col in enumerate(columns):
        p=float(2*t.sf(abs(b[i]/se[i]),G-1)) if se[i]>0 else None
        terms[col]={"beta":float(b[i]),"island_cluster_SE":float(se[i]),
             "nominal_95CI":[float(b[i]-cutoff*se[i]),float(b[i]+cutoff*se[i])],
             "nominal_unadjusted_two_sided_p":p}
    return {
       "status":"estimated_exploratory_photo_source_composition_not_pollen_service",
       "source_GSHHG_independent_islands":int(G),
       "source_island_year_FEs":int(H),
       "true_recorded_photo_grid_years":n,
       "y_recorded_response":y,"weighting":weighting,"predictors":columns,
       "source_predictor_condition_number":cond,
       "original_island_clustered_95CI_nominal":True,
       "terms":terms,"no_insect_moth_abundance_or_effective_pollination":True,
       "observer_effort_and_land_cover_not_fully_measured":True,
    }

def run(source_csv,source_receipt,patch_csv,patch_receipt,out):
    x,rec=eligible_full(source_csv,source_receipt,patch_csv,patch_receipt)
    result={}
    for threshold in THRESHOLDS:
        selected,support=source_stable_group_selection(x,threshold,.05,False)
        selected["centre_lat_c"]=selected.centre_lat-selected.centre_lat.median()
        y="recorded_moth_family_share_pp"
        model=fit_source(selected,y,"island_equal",list(REGISTRY))
        base=fit_source(selected,y,"island_equal",
                        ["prior_log1p_PCNL","log_all_fourfamily_observations","centre_lat_c"])
        photo_control=fit_source(selected,"log_all_fourfamily_observations",
                        "island_equal",["prior_log1p_PCNL","dark_land_fraction",
                                        "largest_dark_component_fraction_of_dark",
                                        "centre_lat_c"])
        source_site_equal=fit_source(selected,y,"site_equal",list(REGISTRY))
        st,sup_strict=source_stable_group_selection(x,threshold,.05,True)
        st["centre_lat_c"]=st.centre_lat-st.centre_lat.median()
        strict=fit_source(st,y,"island_equal",list(REGISTRY))
        result[str(threshold)]={
           "source_eligibility":support,
           "primary_connected_dark_given_light_and_dark_fraction":model,
           "same_exact_source_sites_intensity_without_dark_geometry":base,
           "negative_control_recorded_insect_photo_effort":photo_control,
           "site_equal_documentation_sensitivity":source_site_equal,
           "require_independent_dark_area_span_0p05":{"source_support":sup_strict,
                   "fit":strict}}
    output={
       "schema":SCHEMA,
       "source_previous_year_PCNL_years":list(range(2012,2025)),
       "original_GBIF_source_date_years":list(range(2013,2026)),
       "original_2026_08_GBIF_4families_DOI":"10.15468/dl.8yrbe7",
       "actual_corrected_PCNl_DOI":"10.5281/zenodo.17013414",
       "all_true_original_62688_source_photo_cells":62688,
       "real_photo_source_islands_before_dark_pixel_support":196,
       "prespecified_exploratory_PCNl_dark_thresholds":list(THRESHOLDS),
       "source_photo_records_not_actual_moth_density":True,
       "no_lux_at_flower_no_effective_moth_pollen_deposition":True,
       "source_geo_island_pixels_5km_masked_sea_and_lake":True,
       "raster_dark_fraction_correlated_with_pixel_light_and_contiguity":True,
       "actual_source_and_dark_buffer_threshold_outputs":result}
    out=Path(out)
    require(not out.exists(),"Do not replace original observational result")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("ALAN_GLOBAL_TRUE_DARK_PATCH_INDEPENDENT_MOTH_RECORD_SOURCE",
          json.dumps(output,sort_keys=True),flush=True)
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for field in ("source","source_receipt","patch","patch_receipt","out"):
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.source_receipt,a.patch,a.patch_receipt,a.out)
