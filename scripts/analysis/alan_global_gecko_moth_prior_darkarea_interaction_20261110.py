#!/usr/bin/env python3
"""Global Gekkonidae-source × genuine prior-year PCNL dark-land moderation.

Combine (1) authenticated GBIF gecko documentation preceding insect photo
year, (2) 196-island original moth-versus-butterfly GBIF *photo-record*
composition, and (3) original 2012–2024 GSHHG land-only 5km PCNL dark
fraction/connectivity measured before record years. Positive photo records
do NOT identify gecko feeding, moth abundance or effective pollination.

Same island×year fixed intercepts, geographic source/site volumes, cluster
by real islands and equal weights per island. Negative control prior
Scincidae (non-gecko lizard) and prior Hemidactylus photo evidence.
Original 5km patches can overlap; 10km source-grid separation diagnostic
refits source-island-corrected models. No p-hacking via threshold selection.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    fit_source, SCHEMA as DARK_FIT_SCHEMA)
from scripts.analysis.alan_global_GBIF_darkpatch_spatial_robustness_20261110 import (
    fixed_hash_10km_thinning)
from scripts.analysis.alan_global_gecko_moth_historical_source_cosupport_20261110 import (
    SCHEMA as GECKO_SOURCE_SCHEMA)
from scripts.analysis.alan_global_trueGBIF_PCNl_13year_darkpatch_20261110 import (
    SCHEMA as DARK_SOURCE_SCHEMA)
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import (
    sha256)

SCHEMA="ALAN_true_global_GECKO_photo_documentation_x_5km_prior_darkland_moth_GBIF_source_v1"
RESPONSE="recorded_moth_family_share_pp"
THRESHOLDS=(0.,3.)
BASE={"prior_log1p_PCNL","dark_land_fraction",
      "largest_dark_component_fraction_of_dark",
      "log_all_fourfamily_observations","centre_lat_c"}
TREATMENTS=("prior_gecko_documented","prior_scink_documented",
            "prior_Hemidactylus_documented")

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GECKO_DARK_SOURCE: "+str(msg))

def joined_source(insect,igr,patch,dgr):
    a=json.loads(Path(igr).read_text())
    b=json.loads(Path(dgr).read_text())
    require(a["schema"]==GECKO_SOURCE_SCHEMA and
            a["original_surveyed_moth_insect_geogrid_year_rows"]==62688 and
            sha256(insect)==a["original_source_joint_SHA256"],
            "Authentic gecko/moth source SHA and 196 island IDs not verified")
    require(b["schema"]==DARK_SOURCE_SCHEMA and
            b["source_actual_global_GBIF_photo_cells"]==62688 and
            b["source_original_islands"]==196 and
            sha256(patch)==b["source_2threshold_global_patch_CSV_SHA256"],
            "Prior 13-year PCNL dark-patch source not verified")
    orig=pd.read_csv(insect,dtype={"island":str})
    d=pd.read_csv(patch,dtype={"island":str})
    k=["island","year","latbin","lonbin"]
    require(len(orig)==62688 and
            not orig.duplicated(k).any() and
            len(d)==2*62688 and
            not d.duplicated(k+["threshold_PCNL_native"]).any() and
            set(d.threshold_PCNL_native)==set(THRESHOLDS),
            "Original GECKO or dark source 5km threshold repeated")
    require(orig.prior_log1p_PCNL.notna().all() and orig.year.between(2013,2025).all(),
            "Source photo-year and prior-year light missing")
    extras=["threshold_PCNL_native","status","n_dark_raster_pixels",
            "dark_land_fraction","largest_dark_component_fraction_of_dark",
            "local_mean_PCNL_native","n_valid_original_PCNL_pixels",
            "valid_land_fraction"]
    out=d[k+extras].merge(orig,on=k,how="left",validate="many_to_one",indicator="_dark_geo_match")
    require(out._dark_geo_match.eq("both").all(),
            "Gecko photo sites not identical to actual year-prior dark mosaic sites")
    for name in TREATMENTS:
        out[name]=out[name].astype(str).str.lower().eq("true")
    return out,a,b

def usable(d,t,years=(2015,2025),count_min=5,
           require_dark_span=.05,require_connection_span=.05):
    z=d.loc[d.threshold_PCNL_native.eq(t)&
            d.year.between(years[0],years[1])&
            d.status.eq("observed_satellite_geometry_only")&
            d.n_dark_raster_pixels.gt(0)&
            d.dark_land_fraction.notna()&
            d.largest_dark_component_fraction_of_dark.notna()&
            d.all_four_family_records.ge(count_min)].copy()
    require(z.prior_log1p_PCNL.notna().all(),"No source light after geodesic patch")
    groups=z.groupby("island_year",observed=True).agg(
        n_geo=("lonbin","size"),
        n_luxproxy=("prior_log1p_PCNL","nunique"),
        light_min=("prior_log1p_PCNL","min"),
        light_max=("prior_log1p_PCNL","max"),
        dark_min=("dark_land_fraction","min"),
        dark_max=("dark_land_fraction","max"),
        connect_min=("largest_dark_component_fraction_of_dark","min"),
        connect_max=("largest_dark_component_fraction_of_dark","max"))
    selected=groups.index[(groups.n_geo>=2)&
         (groups.n_luxproxy>=2)&(groups.light_max-groups.light_min>=.02)&
         (groups.dark_max-groups.dark_min>=require_dark_span)&
         (groups.connect_max-groups.connect_min>=require_connection_span)]
    z=z.loc[z.island_year.isin(selected)].copy()
    stats={
        "source_threshold_PCNL_native":t,
        "source_observation_year_window":list(years),
        "source_min_actual_same_cell_insect_records":count_min,
        "source_min_dark_land_fraction_span_per_islandyear":require_dark_span,
        "source_min_dark_connectedness_span_per_islandyear":require_connection_span,
        "true_eligible_geographic_photo_cells":len(z),
        "true_eligible_independent_GSHHG_islands":int(z.island.nunique()),
        "eligible_island_year_fixed_effect_groups":int(z.island_year.nunique()),
        "recorded_prior_gecko_positive_source_islands":int(
            z.loc[z.prior_gecko_documented,"island"].nunique()),
        "recorded_prior_gecko_undocumented_source_islands":int(
            z.loc[~z.prior_gecko_documented,"island"].nunique()),
        "recorded_gecko_source_not_confirmed_biological_occupancy":True,
        "no_original_nighttime_behaviour_or_pollen_transfer":True}
    return z,stats

def source_regression(x,moderator,weight="island_equal",
                      light_gecko_interaction=True,
                      dark_gecko_interaction=True):
    Gpos=int(x.loc[x[moderator],"island"].nunique())
    Gneg=int(x.loc[~x[moderator],"island"].nunique())
    Hpos=int(x.loc[x[moderator],"island_year"].nunique())
    Hneg=int(x.loc[~x[moderator],"island_year"].nunique())
    if Gpos<20 or Gneg<20 or Hpos<40 or Hneg<40:
        return {"status":"not_estimable_unbalanced_prior_reptile_photo_source",
            "genuine_islands_with_prior_gecko_source":Gpos,
            "genuine_islands_without_prior_gecko_photo_source":Gneg,
            "documented_and_undocumented_groups":[Hpos,Hneg],
            "no_gecko_absence_inference":True}
    x=x.copy()
    x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
    cols=["prior_log1p_PCNL","dark_land_fraction",
          "largest_dark_component_fraction_of_dark",
          "log_all_fourfamily_observations","centre_lat_c"]
    if light_gecko_interaction:
        x["moderator_x_prior_light"]=x[moderator].astype(float)*x.prior_log1p_PCNL
        cols.append("moderator_x_prior_light")
    if dark_gecko_interaction:
        x["moderator_x_dark_land"]=x[moderator].astype(float)*x.dark_land_fraction
        cols.append("moderator_x_dark_land")
    r=fit_source(x,RESPONSE,weight,cols)
    return {"status":r.get("status"),"moderator":moderator,
            "real_documented_gecko_source_islands":Gpos,
            "true_undocumented_gecko_photo_islands":Gneg,
            "true_positive_prior_photo_source_island_years":Hpos,
            "true_no_prior_photo_source_island_years":Hneg,
            "original_lux_and_fragmentation_become_modifiers_NOT_predation":True,
            "original_source_group_fixed_intercepts":True,
            "gecko_photo_source_not_true_gecko_absence_or_occupancy":True,
            "source_GBIF_photo_response_not_effective_moth_pollination":True,
            "fit":r}

def run(gecko,gecko_receipt,dark,dark_receipt,out):
    d,a,b=joined_source(gecko,gecko_receipt,dark,dark_receipt)
    reports={}
    for t in THRESHOLDS:
        x,eligibility=usable(d,t)
        measures={
            "main_exploratory_gecko_documentation_x_dark_and_x_light":
                source_regression(x,"prior_gecko_documented"),
            "independent_moth_source_control_no_gecko_interaction":
                fit_source(x.assign(
                   centre_lat_c=x.centre_lat-x.centre_lat.median()),
                   RESPONSE,"island_equal",
                   ["prior_log1p_PCNL","dark_land_fraction",
                   "largest_dark_component_fraction_of_dark",
                   "log_all_fourfamily_observations","centre_lat_c"]),
            "source_Scincidae_photo_prior_comparator":
                source_regression(x,"prior_scink_documented"),
            "source_Hemidactylus_prior_specific":
                source_regression(x,"prior_Hemidactylus_documented"),
            "gecko_x_dark_only_no_light_interaction":
                source_regression(x,"prior_gecko_documented",
                                  light_gecko_interaction=False),
            "island_year_photo_cell_equal_weight":
                source_regression(x,"prior_gecko_documented",weight="site_equal"),
            "source_2019plus":
                source_regression(usable(x,t,years=(2019,2025))[0],
                                  "prior_gecko_documented")
        }
        if t==0.:
            separated={}
            for salt in ("sourcegeo0","sourcegeo1"):
                thinned=fixed_hash_10km_thinning(x,salt=salt,min_km=10.)
                thinned,support=usable(thinned,t)
                separated[salt]={"source_support":support,
                                 "refit":source_regression(thinned,"prior_gecko_documented")}
            measures["source_10km_nonoverlap_5km_geodesic_buffers"]=separated
        reports[str(t)]={"source_eligibility":eligibility,"models":measures}
        print("ALAN_GLOBAL_TRUE_GECKO_MOTH_DARK_AREA_MODERATION",json.dumps({
            "threshold":t,"eligible":eligibility,
            "main":measures["main_exploratory_gecko_documentation_x_dark_and_x_light"],
            "skink":measures["source_Scincidae_photo_prior_comparator"]},
            sort_keys=True),flush=True)
    report={
       "schema":SCHEMA,"actual_original_photo_site_year_source":62688,
       "source_GBIF_reptile_and_moth_snapshot_DOI":"10.15468/dl.8yrbe7",
       "true_preobserved_13year_PCNL_DOI":"10.5281/zenodo.17013414",
       "original_gecko_moth_join_SHA256":sha256(gecko),
       "original_13year_darkpatch_SRGIS_SHA256":sha256(dark),
       "actual_source_contrasts":reports,
       "never_gecko_absence_from_no_GBIF_photo":True,
       "no_natural_moth_population_or_eaten_moth_diet_measured":True,
       "no_effective_night_pollen_transfer_or_plant_seed_set":True,
       "all_pvalues_exploratory_and_multiplicity_unadjusted":True,
       "native_versus_invasive_gecko_status_not_verified_by_island":True}
    out=Path(out);require(not out.exists(),"Do not replace authentic global predator source model")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_GLOBAL_ACTUAL_GECKO_MOTH_DARK_AREA_SOURCE_MODEL_COMPLETE",
          json.dumps({"all_thresholds":list(reports),
                       "original_photo_sites":62688},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for field in ("gecko","gecko_receipt","dark","dark_receipt","out"):
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.gecko,a.gecko_receipt,a.dark,a.dark_receipt,a.out)
