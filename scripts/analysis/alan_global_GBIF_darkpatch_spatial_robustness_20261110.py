#!/usr/bin/env python3
"""Global 2026-08 GBIF moth-PHOTO-source nightlight/dark-area robustness.

Real source 196 island / 62688 photo cells and actual 2012-24 annual PCNL
geodesic 5km source. Sensitivities intentionally reject nearby overlapping
5km circles by choosing 10km-separated source photo-cell centres within each
original island-year, with SHA-sorted selection unrelated to moth outcome,
radiance, dark area or connectivity. All results remain PHOTO DOCUMENTATION,
not real insect population density, flower visits or pollen deposition.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    eligible_full,source_stable_group_selection,fit_source,SCHEMA as CORE_SCHEMA,
    THRESHOLDS,REGISTRY,require)

SCHEMA="ALAN_2026_11_true_global_moth_dark_fractions_spatial_cluster_bias_robustness_v1"
R_EARTH_M=6371008.8

def core_term(result,term="dark_land_fraction"):
    val=result.get("terms",{}).get(term,{})
    return {"status":result.get("status"),
            "n_true_islands":result.get("source_GSHHG_independent_islands",
                    result.get("source_islands")),
            "n_original_photo_grid_years":result.get("true_recorded_photo_grid_years",
                    result.get("source_siteyear_cells")),
            "n_island_year_groups":result.get("source_island_year_FEs",
                    result.get("source_year_groups")),
            "coefficient":val.get("beta"),
            "original_island_cluster_95CI":val.get("nominal_95CI"),
            "nominal_unadjusted_two_sided_p":val.get("nominal_unadjusted_two_sided_p"),
            "mean_light_coefficient":result.get("terms",{}).get(
                    "prior_log1p_PCNL",{}).get("beta"),
            "connected_dark_coefficient":result.get("terms",{}).get(
                    "largest_dark_component_fraction_of_dark",{}).get("beta")}

def ecef(lat,lon):
    rad_lat=np.deg2rad(np.asarray(lat,dtype=float))
    rad_lon=np.deg2rad(np.asarray(lon,dtype=float))
    xyz=np.column_stack([np.cos(rad_lat)*np.cos(rad_lon),
                         np.cos(rad_lat)*np.sin(rad_lon),
                         np.sin(rad_lat)])*R_EARTH_M
    return xyz

def fixed_hash_10km_thinning(df,salt="sourceblinded0",min_km=10.):
    """Greedy geodesic-chord independent set, ordering fixed by original IDs.

    This is NOT a random population sample or new biological replicate:
    use all three declared salts as reproducibility sensitivities. Original
    source site/year identity is retained, and selections use no outcomes.
    """
    require(min_km>=10. and min_km<=40.,"Dangerously overlapping 5km exposure buffers")
    q=df.copy()
    require(not q.duplicated(["island","year","latbin","lonbin"]).any(),
            "Original site/year repeated: impossible geodesic thinning")
    center_lon=.05*(q.lonbin.to_numpy(float)+.5)-180.
    q["_source_lon"]=center_lon
    source_ids=q.index.to_numpy()
    selected_indices=[]
    chord=2*R_EARTH_M*np.sin((min_km*1000)/(2*R_EARTH_M))
    for _,sub in q.groupby(["island","year"],sort=True):
        locations=ecef(sub.centre_lat.to_numpy(float),
                       sub._source_lon.to_numpy(float))
        ids=sub.index.to_numpy()
        hashes=[hashlib.sha256(
                 f"{salt}|{v.island}|{v.year}|{v.latbin}|{v.lonbin}".encode()
                ).hexdigest()
                for v in sub.itertuples(index=False)]
        order=sorted(range(len(sub)),key=lambda i:hashes[i])
        tree=cKDTree(locations)
        usable=np.ones(len(sub),dtype=bool)
        for pos in order:
            if not usable[pos]:continue
            selected_indices.append(ids[pos])
            neighbours=tree.query_ball_point(locations[pos],r=chord)
            usable[neighbours]=False
    out=df.loc[selected_indices].copy()
    require(not out.duplicated(["island","year","latbin","lonbin"]).any(),
            "New geographic source thinning introduced duplicate photo cells")
    return out

def true_fit(x,threshold):
    z,support=source_stable_group_selection(x,threshold,.05,False)
    if len(z):
        z["centre_lat_c"]=z.centre_lat-z.centre_lat.median()
    terms=fit_source(z,"recorded_moth_family_share_pp",
                     "island_equal",list(REGISTRY))
    return {"true_spatial_eligibility":support,
            "source_record_darkland_association":core_term(terms),
            "conditional_connectivity":core_term(
                        terms,"largest_dark_component_fraction_of_dark"),
            "source_mean_light":core_term(terms,"prior_log1p_PCNL")}

def run(source,source_receipt,patch,patch_receipt,model_json,out):
    model=json.loads(Path(model_json).read_text())
    require(model.get("schema")==CORE_SCHEMA and
            model.get("all_true_original_62688_source_photo_cells")==62688 and
            set(model.get("actual_source_and_dark_buffer_threshold_outputs",{}))
            =={"0.0","3.0"},
            "Original jointly sourced PCNL patch/photo source fit not checked")
    all_source,_=eligible_full(source,source_receipt,patch,patch_receipt)
    source_results={}
    for threshold in THRESHOLDS:
        dd=all_source.loc[all_source.threshold_PCNL_native.eq(threshold)].copy()
        require(len(dd)==62688,"True original yearly photo-cell cohort incomplete")
        main=dd.loc[dd.status.eq("observed_satellite_geometry_only") &
                    dd.n_dark_raster_pixels.gt(0) &
                    dd.largest_dark_component_fraction_of_dark.notna()].copy()
        result={"original_pixel_supported_and_real_dark_source_cells":len(main),
                "reprocessed_full_original_main":true_fit(dd,threshold),
                "source_year_sensitivities":{},
                "observational_source_effort_20plus":true_fit(
                  dd.loc[dd.original_total_fourfamily_source_records.ge(20)],
                  threshold),
                "exclude_five_most_photo_dense_islands":[],
                "nonoverlapping_5km_radius_photo_grid_centers_10km_apart":{}}
        for y in (2017,2020,2023):
            result["source_year_sensitivities"][str(y)]=true_fit(
                 dd.loc[dd.year.ge(y)],threshold)
        most=main.groupby("island").size().sort_values(ascending=False).head(5)
        for island,n in most.items():
            sub=dd.loc[~dd.island.eq(island)].copy()
            result["exclude_five_most_photo_dense_islands"].append({
               "omitted_island_GSHHG_source_id":str(island),
               "omitted_documented_photo_grid_years":int(n),
               "re_fitted":true_fit(sub,threshold)})
        for salt in ("sourceblinded0","sourceblinded1","sourceblinded2"):
            # Select only true spatially measured source centres. Originally
            # recorded moth/cell counts are NOT used in selection or sorting.
            sparse=fixed_hash_10km_thinning(main,salt,10.)
            r=true_fit(sparse,threshold)
            result["nonoverlapping_5km_radius_photo_grid_centers_10km_apart"][salt]={
                "one_original_source_site_islandyear_per_10km_proximity_set":True,
                "original_photo_cells_retained_before_regression_support":len(sparse),
                "retained_source_islands_before_regression":int(sparse.island.nunique()),
                "re_fitted":r}
        source_results[str(threshold)]=result
        print("ALAN_SOURCE_TRUE_GLOBAL_DARK_AREA_GEOGRAPHIC_STRESS",json.dumps({
           "threshold":threshold,
           "main":result["reprocessed_full_original_main"]["source_record_darkland_association"],
           "spatial_separation_10km":{salt:obj["re_fitted"]["source_record_darkland_association"]
                                     for salt,obj in result[
                                      "nonoverlapping_5km_radius_photo_grid_centers_10km_apart"].items()},
           "2017_2020_2023":{k:v["source_record_darkland_association"] for k,v in
                               result["source_year_sensitivities"].items()}},
           sort_keys=True),flush=True)
    report={
       "schema":SCHEMA,"source_moth_GBIF_photo_record_NOT_abundance_or_pollination":True,
       "original_true_2026_08_GBIF_DOI":"10.15468/dl.8yrbe7",
       "original_corrected_PCNL_prior_year_DOI":"10.5281/zenodo.17013414",
       "source_actual_island_group_fixed_effect_and_original_island_clusters":True,
       "source_genuine_original_2013to2025_photo_cells_prior_dark_62688":True,
       "spatial_thinning_never_based_on_moth_outcome_ALAN_or_dark_connectivity":True,
       "original_5km_buffer_overlap_removed_with_10km_grid_point_spacing":True,
       "source_islands_with_source_missing_month_photographer_or_local_lux":True,
       "source_negative_photo_sampling_control_not_measured_moth_population":True,
       "all_output_coefficient_claims_exploratory_not_confirmatory":True,
       "source_threshold_and_geographic_robustness":source_results}
    out=Path(out)
    require(not out.exists(),"No replacement of original signed sensitivity")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_GLOBAL_DARK_AREA_SPATIAL_ROBUSTNESS_COMPLETE",
          json.dumps({"source_valid":True,
                      "thresholds":list(source_results)},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for field in ("source","source_receipt","patch","patch_receipt","model_json","out"):
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.source_receipt,a.patch,a.patch_receipt,a.model_json,a.out)
