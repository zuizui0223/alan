#!/usr/bin/env python3
"""Gecko-vs-skink SOURCE photographed proportions under real solar day vs night.

Original GBIF 2026-08 eventtime with explicit UTC offset (not generic
date-only eventDate), solar altitude at original point, original GSHHG
island location. Record groups compare gecko v Scincidae photographic
observations in same recorded SOLAR class and same island×observation
year. Previous-year PCNL 2012-24 matched to exactly same ~0.05° cell
centers. Not actual time-integrated gecko activity, trophic predation,
moth abundance or effective nocturnal pollen services.

If same-island solar-day or solar-night reptile photo classes lack at
least 30 independent islands and 2 source points with real light contrasts,
report original support rather than fitting a spurious slope.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_gecko_skink_recorded_light_niche_FE_20261110 import (
    load_original,require)
from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import fit_source
from scripts.analysis.alan_true_GECKO_GBIF_solar_daynight_source_20261110 import (
    GECKO_FAMILIES,SCHEMA as CLOCK_SCHEMA,sha256)

SCHEMA="ALAN_GBIF_gecko_skink_solar_day_vs_night_prior_PCNL_record_niche_v1"
CLASSES=("solar_night","solar_day")
Y="recorded_gecko_photo_fraction_pp"

def prepare(true_time,time_receipt,true_light,light_receipt):
    r=json.loads(Path(time_receipt).read_text())
    require(r.get("schema")==CLOCK_SCHEMA and
            r.get("original_full_Squamata_source_rows")==11634053 and
            sha256(true_time)==r["source_solar_photo_CSV_SHA256"],
            "Not signed genuine gecko solar night/day source")
    time=pd.read_csv(true_time,dtype={"island":str},low_memory=False)
    keys=["island","year","latbin","lonbin"]
    need=set(keys+["source_taxon_family","verified_solar_photo_class",
                   "source_GBIF_photo_record_count"])
    require(need.issubset(time.columns) and time.year.between(2013,2025).all() and
            set(time.source_taxon_family).issubset(set(GECKO_FAMILIES)|{"Scincidae"}) and
            not time.duplicated(keys+["source_taxon_family","verified_solar_photo_class"]).any(),
            "Solar time original photo grid source keys duplicated")
    time["source_gecko_photo"]=np.where(time.source_taxon_family.isin(GECKO_FAMILIES),
                                   time.source_GBIF_photo_record_count,0)
    time["source_scink_photo"]=np.where(time.source_taxon_family.eq("Scincidae"),
                                   time.source_GBIF_photo_record_count,0)
    d=time.groupby(keys+["verified_solar_photo_class"],sort=True)[
        ["source_gecko_photo","source_scink_photo"]].sum().reset_index()
    require(int(d.source_gecko_photo.sum())==
            r["all_gecko_family_island_observations_with_valid_solar_clock"] and
            int(d.source_scink_photo.sum())==
            r["all_scink_island_observations_with_valid_solar_clock"],
            "True source solar time island photo totals drifted")
    source,_=load_original(true_light,light_receipt)
    require(not source.duplicated(keys).any() and len(source)==34056,
            "Old 2026-08 source gecko site-year recording frame altered")
    joined=d.merge(source[keys+["log1p_prior_PCNL","latitude_center",
                         "area_km2","original_total_reptile_photo_source_events"]],
                   on=keys,how="left",validate="many_to_one",indicator="_solar_gbif_join"
    require(joined._solar_gbif_join.eq("both").all() and
            joined.source_gecko_photo.sum()>0 and
            joined.source_scink_photo.sum()>0,
            "Actual timestamped GBIF source is not a subset of original geolocated reptile photos")
    joined["source_total_solar_class_photos"]=(
        joined.source_gecko_photo+joined.source_scink_photo)
    joined[Y]=100*joined.source_gecko_photo/joined.source_total_solar_class_photos
    joined["source_log_photo_events"]=np.log1p(joined.source_total_solar_class_photos)
    joined["centre_lat"]=joined.latitude_center
    return joined,r

def cell_eligibility(data,solar_class,minimum_each=1):
    z=data.loc[data.verified_solar_photo_class.eq(solar_class)&
               data.log1p_prior_PCNL.notna()].copy()
    z["island_year"]=z.island.astype(str)+"|"+z.year.astype(str)
    gs=z.groupby("island_year",sort=True).agg(
        island=("island","first"),n_photo_cells=("lonbin","size"),
        n_light_levels=("log1p_prior_PCNL","nunique"),
        min_light=("log1p_prior_PCNL","min"),
        max_light=("log1p_prior_PCNL","max"),
        gecko_photos=("source_gecko_photo","sum"),
        scink_photos=("source_scink_photo","sum"))
    valid=gs.index[(gs.n_photo_cells>=2)&(gs.n_light_levels>=2)&
                   (gs.max_light-gs.min_light>=.02)&
                   gs.gecko_photos.ge(minimum_each)&
                   gs.scink_photos.ge(minimum_each)]
    z=z.loc[z.island_year.isin(valid)].copy()
    stat={
       "solar_class":solar_class,"min_source_photo_records_either_reptile_class":minimum_each,
       "original_exact_timezone_offset_source_cells_with_numeric_preyear_light":int(
           data.loc[data.verified_solar_photo_class.eq(solar_class),
                    "log1p_prior_PCNL"].notna().sum()),
       "supported_island_year_groups":len(valid),
       "supported_true_GSHHG_islands":int(z.island.nunique()),
       "supported_true_reptile_observation_source_cells":len(z),
       "actual_source_gecko_photos":int(z.source_gecko_photo.sum()),
       "actual_source_scink_photos":int(z.source_scink_photo.sum()),
       "recorded_no_photo_is_not_real_reptile_absence":True}
    return z,stat

def regression(z,weight):
    if not len(z):return {"status":"not_estimable_no_solar_class_photo_coobservations"}
    z=z.copy()
    z["prior_log1p_PCNL"]=z.log1p_prior_PCNL
    z["centre_lat_c"]=z.centre_lat-z.centre_lat.median()
    z["log_all_fourfamily_observations"]=z.source_log_photo_events
    return fit_source(z,Y,weight,
        ["prior_log1p_PCNL","log_all_fourfamily_observations","centre_lat_c"])

def run(timed,time_receipt,light,light_receipt,out):
    original,provenance=prepare(timed,time_receipt,light,light_receipt)
    comparisons={}
    for category in CLASSES:
        selected,support=cell_eligibility(original,category,1)
        strict,support_strict=cell_eligibility(original,category,5)
        comparisons[category]={
          "original_source_eligibility":support,
          "island_equal_documented_reptile_families":regression(selected,"island_equal"),
          "source_cell_year_equal_sensitivity":regression(selected,"site_equal"),
          "at_least_five_gecko_and_skink_photos_each_island_year":{
              "source_support":support_strict,
              "fit":regression(strict,"island_equal")}}
        print("ALAN_ACTUAL_ISLAND_GECKO_SUN_OBSERVED_PHOTO_LIGHT",
            json.dumps({"solar_class":category,
              "original_support":support,
              "primary":comparisons[category]["island_equal_documented_reptile_families"]},
              sort_keys=True),flush=True)
    output={
        "schema":SCHEMA,
        "original_GBIF_2026_08_DOI":"10.15468/dl.8yrbe7",
        "source_pcNl_2012to24_DOI":"10.5281/zenodo.17013414",
        "true_source_timed_gecko_scink_GSHHG_photo_events":{
          "gecko":provenance["all_gecko_family_island_observations_with_valid_solar_clock"],
          "scink":provenance["all_scink_island_observations_with_valid_solar_clock"]},
        "true_solar_class_source_photo_comparisons":comparisons,
        "no_true_gecko_moth_predation_or_moth_pollen_transfer":True,
        "observed_photo_clock_not_animal_time_budget":True,
        "clock_timezone_present_but_may_have_provider_errors":True,
        "all_source_models_exploratory_not_confirmatory":True}
    out=Path(out);require(not out.exists(),"Cannot replace true gecko sunlight field")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_GECKO_SOLAR_NIGHT_AND_DAY_SOURCE_LIGHT_CHECK_COMPLETE",json.dumps({
         "night":comparisons["solar_night"]["original_source_eligibility"],
         "day":comparisons["solar_day"]["original_source_eligibility"]},
        sort_keys=True),flush=True)
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("timed","time_receipt","light","light_receipt","out"):
        p.add_argument("--"+x.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.timed,a.time_receipt,a.light,a.light_receipt,a.out)
