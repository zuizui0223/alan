#!/usr/bin/env python3
"""GECKO-v-SCINK GBIF PHOTO-RECORD light-niche contrast at actually surveyed sites.

Original 2026-08 GBIF Squamata source: 66,038 gecko family photos, 35,361
Scincidae photos from original 2013–25 dates, same GSHHG>=20km2 island
geographic 0.05deg grid/year. Pixel intensity is the real previous calendar
year PCNL 2012–24; missing/lake/sea-centered source cells are excluded.

Key design improvement: REQUIRE BOTH gecko and skink photo categories
documented within SAME island × observation year to avoid "gecko absent"
labels. Still not a standardized lizard census and photographic activity
bias remains. Response is SOURCE DOCUMENTATION percentage of gecko family
among gecko+skink human observation records, never true occupancy, gecko
predation on moths, nocturnal activity or vegetation pollen service.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_GBIF_moth_record_fraction_within_island_FE_20261110 import (
    original_group_demean,transform_weights)
from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    fit_source)
from scripts.analysis.alan_global_GBIF_darkpatch_spatial_robustness_20261110 import (
    fixed_hash_10km_thinning)
from scripts.analysis.alan_gecko_skink_original_13year_priorPCNL_20261110 import (
    SCHEMA as SOURCE_SCHEMA,sha256,require)

SCHEMA="ALAN_2026_11_gecko_v_scink_GBIF_photo_records_prior_PCNL_same_islandyear_v1"
RESPONSE="gecko_photo_record_percent"
RADIANCE="log1p_prior_PCNL"
PER_GROUP_MIN_LIGHT_SPAN=.02
MIN_GROUP_TAXON_PHOTOS=1

def load_original(csv,receipt):
    d=json.loads(Path(receipt).read_text())
    require(d.get("schema")==SOURCE_SCHEMA and
            d.get("all_original_gecko_and_scink_photo_record_events")==101399 and
            d.get("original_Gekkonidae_plus_other5photo_event_count")==66038 and
            d.get("original_Scincidae_photo_event_count")==35361 and
            sha256(csv)==d.get("original_reptile_photo_cell_CSV_SHA256"),
            "Not actual original 13-year gecko+skink photo raster source")
    z=pd.read_csv(csv,dtype={"island":str},low_memory=False)
    k=["island","year","latbin","lonbin"]
    require(not z.duplicated(k).any() and z.year.between(2013,2025).all() and
            int(z.gecko_source_photo_count.sum())==66038 and
            int(z.scink_source_photo_count.sum())==35361 and
            (z.prior_PCNL_year==z.year-1).all(),
            "Original 101399 reptile family photo source cannot reconcile site/year")
    z[RESPONSE]=100*z.gecko_source_photo_count/z.original_total_reptile_photo_source_events
    z["log_original_reptile_record_count"]=np.log1p(
        z.original_total_reptile_photo_source_events)
    z["centre_lat"]=z.latitude_center
    z["island_year"]=z.island.astype(str)+"|"+z.year.astype(str)
    require(z[RESPONSE].between(0,100).all(),"Source photo fraction impossible")
    return z,d

def support_coobserved(d,min_each=1,min_cell_photo_events=1,
                       earliest_year=2013,
                       min_sameyear_light_span=PER_GROUP_MIN_LIGHT_SPAN):
    require(min_each>=1 and min_cell_photo_events>=1,"Need real source record category")
    z=d.loc[d.year.ge(earliest_year)&
            d.log1p_prior_PCNL.notna()&
            d.original_total_reptile_photo_source_events.ge(min_cell_photo_events)].copy()
    groups=z.groupby("island_year",sort=True).agg(
        island=("island","first"),
        actually_observed_photo_cells=("lonbin","size"),
        source_distinct_light_values=(RADIANCE,"nunique"),
        source_min_light=(RADIANCE,"min"),
        source_max_light=(RADIANCE,"max"),
        source_gecko_records=("gecko_source_photo_count","sum"),
        source_skink_records=("scink_source_photo_count","sum"))
    valid=groups.index[
        groups.actually_observed_photo_cells.ge(2)&
        groups.source_distinct_light_values.ge(2)&
        (groups.source_max_light-groups.source_min_light).ge(min_sameyear_light_span)&
        groups.source_gecko_records.ge(min_each)&
        groups.source_skink_records.ge(min_each)]
    z=z.loc[z.island_year.isin(valid)].copy()
    g=groups.loc[valid]
    summary={
        "actual_within_islandyear_BOTH_photo_categories_observed":True,
        "min_actual_photo_record_per_family_in_year":min_each,
        "min_actual_photo_source_events_per_cell":min_cell_photo_events,
        "source_photos_earliest_observed_year":earliest_year,
        "original_yearminus1_light_min_span":min_sameyear_light_span,
        "original_source_0p05_reptile_cellyears_before_filter":len(d),
        "source_true_PCNL_numeric_reptile_cellyears":int(
            d.log1p_prior_PCNL.notna().sum()),
        "source_retained_grid_cell_years":len(z),
        "source_retained_independent_islands":int(z.island.nunique()),
        "source_retained_observed_island_year_groups":len(g),
        "source_retained_gecko_photo_events":int(z.gecko_source_photo_count.sum()),
        "source_retained_skink_photo_events":int(z.scink_source_photo_count.sum()),
        "no_photo_record_in_cell_means_NOT_DOCUMENTED_not_biological_absence":True}
    return z,summary

def true_wls(d,weight="island_equal",include_richness_proxy=True,
             include_area_slope=False,allowed_island=None):
    x=d.copy()
    if allowed_island is not None:
        x=x.loc[x.island.ne(allowed_island)].copy()
    if x.empty:
        return {"status":"not_estimable_empty_authenticated_source"}
    x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
    x["log_all_fourfamily_observations"]=x.log_original_reptile_record_count
    x["prior_log1p_PCNL"]=x.log1p_prior_PCNL
    regressors=["prior_log1p_PCNL","centre_lat_c"]
    if include_richness_proxy:regressors.append("log_all_fourfamily_observations")
    if include_area_slope:
        x["island_area_c"]=np.log10(x.area_km2)-np.log10(x.area_km2).median()
        x["island_area_x_source_light"]=x.island_area_c*x.prior_log1p_PCNL
        regressors.append("island_area_x_source_light")
    return fit_source(x,RESPONSE,weight,regressors)

def subset_result(d,min_each=1,min_cell=1,year_min=2013,
                  mode="island_equal",exclude=None):
    z,s=support_coobserved(
        d,min_each=min_each,min_cell_photo_events=min_cell,
        earliest_year=year_min)
    result=true_wls(z,mode,True,False,exclude)
    return {"original_source_eligibility":s,"source_regression":result}

def run(original,receipt,out):
    x,meta=load_original(original,receipt)
    z,eligible=support_coobserved(x)
    most_popular=z.groupby("island").size().sort_values(ascending=False)
    largest=most_popular.index[0] if len(most_popular) else None
    trials={
       "primary_both_reptile_photo_families_same_islandyear":
             subset_result(x),
       "without_log_total_reptile_photo_count":{
            "original_source_eligibility":eligible,
            "source_regression":true_wls(z,"island_equal",False)},
       "siteyear_equal_source_weighting":
             subset_result(x,mode="site_equal"),
       "islandyear_equal_source_weighting":
             subset_result(x,mode="island_year_equal"),
       "same_source_BOTH_families_5plus_photo_events_each_islandyear":
             subset_result(x,min_each=5),
       "source_cell_at_least_3_reptile_photo_records":
             subset_result(x,min_cell=3),
       "source_2017_plus":
             subset_result(x,year_min=2017),
       "source_2020_plus":
             subset_result(x,year_min=2020),
       "source_island_area_x_light_sensitivity":{
           "original_source_eligibility":eligible,
           "source_regression":true_wls(z,"island_equal",True,True)},
       "source_largest_reptile_photo_island_omitted":
             subset_result(x,exclude=largest),
       "negative_control_total_gecko_plus_skink_photo_event_volume":{
           "original_source_eligibility":eligible,
           "source_regression":fit_source(
               z.assign(prior_log1p_PCNL=z.log1p_prior_PCNL,
                        centre_lat_c=z.centre_lat-z.centre_lat.median()),
               "log_original_reptile_record_count","island_equal",
               ["prior_log1p_PCNL","centre_lat_c"])},
       "independently_10km_spaced_source_photo_cells":{}
    }
    for salt in ("reptilegeo0","reptilegeo1"):
        thin=fixed_hash_10km_thinning(z,salt=salt,min_km=10.)
        v,s=support_coobserved(thin)
        trials["independently_10km_spaced_source_photo_cells"][salt]={
            "source_grid_cell_years_pre_contrast":len(thin),
            "original_source_eligibility":s,
            "source_regression":true_wls(v)}
    report={
       "schema":SCHEMA,
       "frozen_2026_08_original_GBIF_snap_DOI":"10.15468/dl.8yrbe7",
       "source_2012to2024_PCNL_original_DOI":"10.5281/zenodo.17013414",
       "real_original_101399_reptile_photo_source_SHA256":sha256(original),
       "main_original_5km_fine_light_resolution_NOT_measured":True,
       "two_source_taxa_compared_both_directly_photo_documented":True,
       "true_gecko_reptile_photo_composition_not_predator_count_or_occupancy":True,
       "all_source_recording_and_detection_biases_unresolved":True,
       "all_source_13year_observational_model_variants":trials,
       "no_streetlight_0to100m_behavioral_measurement":True,
       "gecko_moth_predation_AND_flower_pollen_services_not_observed":True,
       "all_significance_exploratory":True}
    out=Path(out);require(not out.exists(),"Do not replace source-original inference")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_GECKO_V_SCINK_HISTORIC_LIGHT_SOURCE_INFERENCE",
       json.dumps({k:v for k,v in trials.items() if k in (
          "primary_both_reptile_photo_families_same_islandyear",
          "same_source_BOTH_families_5plus_photo_events_each_islandyear",
          "source_largest_reptile_photo_island_omitted",
          "negative_control_total_gecko_plus_skink_photo_event_volume",
          "independently_10km_spaced_source_photo_cells")},
       sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.receipt,a.out)
