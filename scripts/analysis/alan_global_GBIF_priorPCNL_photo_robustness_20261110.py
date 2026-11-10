#!/usr/bin/env python3
"""Postexploratory confidence and source-visibility stress checks for global ALAN.

Uses only already signed 22.7m original GBIF record source, 228361 cell-years
and true PCNL prior-year exposure. Does not fabricate moth/insect counts nor
perform significance-driven cohort choice. Model remains an observation-only
source reporting association, not real insect abundance/efficiency.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_GBIF_moth_record_fraction_within_island_FE_20261110 import (
    source,fit,SCHEMA as MAIN_SCHEMA)
from scripts.analysis.alan_global_GBIF_true4family_priorPCNL_20261110 import (
    sha256,require)
SCHEMA="ALAN_worldwide_GBIF_photo_family_prior_PCNL_source_robustness_v1"
Y="recorded_moth_family_share_pp"

def short(s):
    term=s.get("predictor_coefficients",{}).get("prior_log1p_PCNL",{})
    return {"status":s.get("status"),
            "actual_islands":s.get("observed_islands",s.get("status_counts",{}).get("independent_islands")),
            "actual_site_years":s.get("original_source_sample_size",s.get("status_counts",{}).get("eligible_site_years")),
            "island_year_groups":s.get("island_year_FEs",s.get("status_counts",{}).get("true_island_year_groups")),
            "beta_moth_photo_fraction_pp":term.get("beta"),
            "nominal_95CI_pp":term.get("nominal_95CI"),
            "nominal_two_sided_p":term.get("nominal_two_sided_unadjusted_p")}

def run(source_csv,source_receipt,main_json,out):
    x,p=source(source_csv,source_receipt)
    main=json.loads(Path(main_json).read_text())
    require(main["schema"]==MAIN_SCHEMA and
            main["source_signed_original_source_CSV_SHA256"]==sha256(source_csv) and
            main["real_global_original_source_record_count"]==22713063 and
            main["real_global_original_site_year_cells"]==228361,
            "Previously reported source model receipt cannot be reconciled")
    all_specs=main["all_prespecified_source_model_sensitivities"]
    formatted={name:short(s) for name,s in all_specs.items()}
    years={}
    for y in (2017,2020,2023):
        filt=x.loc[x.year.ge(y)].copy()
        years[str(y)]={"island_equal_source_total_adjusted":
                    short(fit(filt,Y,"island_equal",True,False,False,5)),
                       "same_site_observer_volume_negative_control":
                    short(fit(filt,"log_all_fourfamily_observations",
                              "island_equal",False,False,False,5))}
    # Source-documented top most active island may dominate despite island equal
    # weights because its many independent photo pixels anchor a very different
    # intensity range. Deletion is prespecified rank-based, NOT p-value chosen.
    eligible=x.loc[x.prior_log1p_PCNL.notna() &
                   x.moth_family_records.gt(0) &
                   x.butterfly_family_records.gt(0) &
                   x.all_four_family_records.ge(5)]
    counts=eligible.groupby("island").agg(
        source_cell_years=("year","size"),
        source_recorded_events=("all_four_family_records","sum")).sort_values(
          ["source_cell_years","source_recorded_events"],ascending=False)
    most_popular=counts.index[:5].tolist()
    top_island=[]
    for island in most_popular:
        z=fit(x.loc[~x.island.eq(island)],Y,"island_equal",True,False,False,5)
        top_island.append({"excluded_popular_GSHHG_island_ID":str(island),
                           "original_popular_island_recorded_cell_years":int(
                                  counts.loc[island,"source_cell_years"]),
                           "after_exclusion":short(z)})
    # The strict center-on-land method especially restricts small/coastal islands.
    # Re-estimate in independently source-defined area strata, without
    # claiming random selection or changing the primary.
    area={}
    for name,mask in (
        ("smaller_GSHHG_islands_20_to_500km2",
         x.area_km2.between(20,500,inclusive="left")),
        ("large_GSHHG_islands_500km2_or_more",
         x.area_km2.ge(500)),
        ("higher_same_site_source_event_opportunity_20plus",
         x.all_four_family_records.ge(20)),
        ("lower_source_event_opportunity_5to19",
         x.all_four_family_records.between(5,19))
    ):
        area[name]=short(fit(x.loc[mask],Y,"island_equal",True,False,False,5))
    observed=list(formatted.items())
    final={
        "schema":SCHEMA,
        "original_GBIF_2026_08_DOI":"10.15468/dl.8yrbe7",
        "original_13year_corrected_PCNL_DOI":"10.5281/zenodo.17013414",
        "original_true_global_source_CSV_SHA256":sha256(source_csv),
        "full_source_has_228361_original_GSHHG_islandyearcells":True,
        "original_23million_source_observations_NOT_trap_individuals":True,
        "original_frozen_primary_and_all_sensitivities":formatted,
        "restricted_observation_years":years,
        "most_source_popular_island_omissions":top_island,
        "source_island_area_and_observation_event_strata":area,
        "all_analyses_are_postexploratory_sensitivity_not_confirmatory":True,
        "no_Moth_GloBI_plant_or_Pollen_function_observations":True,
        "no_independent_photographer_or_dataset_detection_effort_normalization":True
    }
    out=Path(out)
    require(not out.exists(),"No overwrite source stress test")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(final,sort_keys=True,indent=2)+"\n")
    print("ALAN_GLOBAL_TRUE_MOTH_SOURCE_LIGHT_ROBUSTNESS",
          json.dumps({"primary":formatted.get("island_equal_primary_recorded_guild_fraction"),
                      "site_equal":formatted.get("site_equal_weighting"),
                      "island_year_equal":formatted.get("island_year_equal_weighting"),
                      "no_effort_adjustment":formatted.get("source_effort_not_adjusted"),
                      "strict_20records":formatted.get("strict_20_source_records_and_both"),
                      "source_zeros_retained":formatted.get(
                           "include_recorded_zero_category_but_not_true_absence"),
                      "years":years,"area_strata":area,
                      "popular_island_deletions":top_island},sort_keys=True),flush=True)
    return final

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--main-model",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.receipt,a.main_model,a.out)
