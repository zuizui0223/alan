#!/usr/bin/env python3
"""Giavi 2021: within-meadow wet-vs-drier floral cover and source selection audit.

Source: Giavi et al Nat Commun 2021 original Zenodo 4540407, original source
RData MD5 f60b21d82c163cbb53520fda0bc61985, exactly 1884 rows,
six physically paired lit and unlit meadows, 37 sampling events.
The published plab proxy measures *surface covered by flowers* along 100 m
transects with ~5-cm circular proxy units. It is NOT a count of flowers,
nor a direct test of ALAN-induced changes in plant reproductive success.

Fresh question: does the apparent negative moisture-vs-lit-dark floral
cover contrast persist when BOTH wet and drier plants are simultaneously
observed in the same meadow and same original sampling event?
Also quantify which plant×site×repeat combinations have only one original
C/T arm, rather than assuming absent plant = zero.

Within each eligible event, average per-plant log1p flower-covered surface
separately for F>7 and F<=7 (within the event, each species equal weight)
then take the wet-minus-drier gap for the lit and dark arm and their
difference. Each independent meadow's event gaps are averaged to report
the number of truly independent lighting contrasts. With only 1-2 meadow
pairs, no population-level wetness moderation is claimed; two-sided
signflip minimum attainable p is 0.5 with 2 meadows.

Also test dry-only F gradient after removing all >7 taxa as falsification
of a purportedly wet-specific light mechanism. F is a species ecological
indicator, NOT directly measured soil or wetland change.
"""
from __future__ import annotations
import argparse, collections, itertools, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_giavi2021_true6meadow_floweradjusted_moisture_20261120 import (
    matched_source, add_moisture, estimate, site_signflip, source_support
)
from scripts.analysis.alan_giavi2021_originalRdata_meadowpair_support_20261120 import fetch

SCHEMA="alan_Giavi2021_same_meadow_same_event_wet_drier_flower_source_gate_v1"
TRAITS=("F_PLANTATT","F_Doring")

def guard(ok,msg):
    if not ok:raise RuntimeError("GIAVI_SAME_MEADOW_WET_DRY_SOURCE_STOP: "+msg)

def event_gaps(data,trait):
    f=data[data[trait].notna()].copy()
    events=[]
    for (site,pair),z in f.groupby(["site","pair_key"],sort=True):
        wet=z[z[trait]>7]
        dry=z[z[trait]<=7]
        if wet.empty or dry.empty:
            continue
        def avg_log(sub,col):
            return float(np.log1p(sub[col].to_numpy(float)).mean())
        lit=avg_log(wet,"flowers_T")-avg_log(dry,"flowers_T")
        dark=avg_log(wet,"flowers_C")-avg_log(dry,"flowers_C")
        source_delta=(lit-dark)
        direct=wet["delta_flower"].mean()-dry["delta_flower"].mean()
        guard(abs(direct-source_delta)<1e-10,
              "event wet-drier original Δflower arithmetic identity violated")
        events.append({
          "original_meadow":site,
          "original_sampling_pair":pair,
          "source_n_wet_plant_taxa":int(wet.plant.nunique()),
          "source_n_drier_plant_taxa":int(dry.plant.nunique()),
          "source_wet_minus_drier_dark_log_flower_cover":dark,
          "source_wet_minus_drier_lit_log_flower_cover":lit,
          "source_double_difference_light_minus_dark_in_wet_drier":source_delta,
          "observed_source_plant_taxa_not_invented":True,
        })
    sites={}
    for site in sorted({x["original_meadow"] for x in events}):
        z=[x for x in events if x["original_meadow"]==site]
        sites[site]={
          "independent_original_meadow_pair":site,
          "sampling_events_with_both_wet_and_drier":len(z),
          "wet_minus_drier_dark_log_cover_event_mean":
              float(np.mean([x["source_wet_minus_drier_dark_log_flower_cover"] for x in z])),
          "wet_minus_drier_lit_log_cover_event_mean":
              float(np.mean([x["source_wet_minus_drier_lit_log_flower_cover"] for x in z])),
          "wet_vs_drier_light_minus_dark_flower_double_difference":
              float(np.mean([x["source_double_difference_light_minus_dark_in_wet_drier"] for x in z]))
        }
    per=[z["wet_vs_drier_light_minus_dark_flower_double_difference"]
         for z in sites.values()]
    return {
      "original_matched_same_site_same_sampling_event_comparisons":len(events),
      "independent_original_meadow_site_pairs_with_qualifying_event":len(sites),
      "independent_meadow_event_averaged_values":sites,
      "source_within_meadow_light_minus_dark_group_difference_site_mean":
          float(np.mean(per)) if per else None,
      "source_direction_negative_meadow_pairs":sum(x<0 for x in per),
      "source_direction_positive_meadow_pairs":sum(x>0 for x in per),
      "exact_two_sided_signflip_min_resolution_if_nonzero_all_sites":
          2/(2**len(per)) if per else None,
      "cannot_infer_wet_group_specific_light_effect":
          len(sites)<3 or len(per)<3,
      "source_event_level_independent_treatment_units_not_assumed":True,
    }

def original_one_arm_selection(obj,names,trait):
    src=obj["DAT"].copy()
    for col in ("site","repeate","plant","treatment"):
        src[col]=src[col].astype(str).str.strip()
    g=src.groupby(["site","repeate","plant"]).treatment.agg(
        lambda x:frozenset(x)).reset_index(name="original_observed_arms")
    guard(len(g)>193,"source arm selection lacks unpaired original observations")
    q=[]
    for z in g.itertuples(index=False):
        info=names.get(z.plant)
        val=info.get(trait) if info else None
        if val is None:continue
        a=z.original_observed_arms
        guard(a in (frozenset({"C"}),frozenset({"T"}),frozenset({"C","T"})),
              "unexpected original source coverage pattern")
        q.append({"site":z.site,"repeate":z.repeate,"plant":z.plant,
                  "moisture":val,"wet":val>7,"arms":a})
    totals={}
    for group,record in [("wet_F_gt7",[x for x in q if x["wet"]]),
                         ("drier_F_le7",[x for x in q if not x["wet"]])]:
        codes=collections.Counter(
            "both" if x["arms"]==frozenset({"C","T"}) else
            "dark_only" if x["arms"]==frozenset({"C"}) else
            "lit_only" for x in record)
        totals[group]={
            "original_source_site_event_plant_keys":len(record),
            "observed_original_lit_and_dark":codes["both"],
            "observed_original_dark_only":codes["dark_only"],
            "observed_original_lit_only":codes["lit_only"],
            "source_both_arm_coverage_fraction":
                codes["both"]/len(record) if record else None}
    return {"original_plant_arm_selection_by_trait":totals,
            "original_unpaired_plant_is_not_assumed_absent_biomass":True,
            "not_habitat_disappearance_or_treatment_extinction_measure":True}

def analyse(source,traits):
    matched,schema=matched_source(source)
    with_trait,names,kew_n=add_moisture(matched)
    result={}
    for ftrait in TRAITS:
        scored=with_trait[with_trait[ftrait].notna()].copy()
        dry=scored[scored[ftrait]<=7].copy()
        wet=scored[scored[ftrait]>7].copy()
        e=event_gaps(scored,ftrait)
        original_arms=original_one_arm_selection(source,names,ftrait)
        dry_only=estimate(dry,ftrait,"flower",False)
        wet_only=estimate(wet,ftrait,"flower",False)
        bysite={}
        for site,z in scored.groupby("site"):
            zz=z[["plant",ftrait]].drop_duplicates()
            # Site-only original within 37 sampled events, each plant
            # must share C/T original source, no artificial zeros.
            slope=estimate(z,ftrait,"flower",False)
            bysite[site]={
                "source_meadow":site,
                "original_matched_scored_plant_pairs":len(z),
                "original_plant_species_scored":len(zz),
                "F_min":float(zz[ftrait].min()),
                "F_max":float(zz[ftrait].max()),
                "wet_species":int(zz[ftrait].gt(7).sum()),
                "drier_species":int(zz[ftrait].le(7).sum()),
                "within_sample_pair_FE_beta_F":
                  slope["beta_F_per_ordinal_point"] if slope else None,
                "inference_independent_of_other_meadows":False
            }
        result[ftrait]={
          "original_source_scored_matched_rows":len(scored),
          "original_wet_plant_names":int(wet.plant.nunique()),
          "original_drier_plant_names":int(dry.plant.nunique()),
          "within_meadow_event_matched_wet_drier_flower":e,
          "original_single_arm_plant_recording_selection":original_arms,
          "dry_only_F_le7_flower_association":dry_only,
          "wet_only_F_gt7_flower_association":wet_only,
          "per_physical_meadow_within_event_source_species_support":bysite,
          "observed_source_FLORAL_PROXY_is_covered_surface_not_flower_number":True,
          "original_light_assignment_not_proven_randomized_by_our_source":True,
          "no_general_wet_specific_moderator_estimated":True
        }
    guard(schema["meadow_pairs"]==6 and schema["matched_plant_pair_cells"]==193,
          "true original Giavi source physical design changed")
    return {
      "schema":SCHEMA,
      "source_original_Giavi_Zenodo_RData_MD5":
          "f60b21d82c163cbb53520fda0bc61985",
      "source_Kew_BIFloraExplorer_Git_blob_SHA1":
          "e6a1da929af81f370d2a1b79177b489e3e573a06",
      "source_original_Giavi_DOI":"10.1038/s41467-021-22011-8",
      "source_original_experimental_meadow_pairs":6,
      "source_sampling_meadow_pair_events":37,
      "source_plant_matched_daytime_flower_visits":193,
      "source_original_plant_and_cover_definition":
          "flower-covered transect surface units of ca 5 cm diameter per circle, measured shortly before visit sampling",
      "source_actual_physical_sites_not_193_independent_treatments":True,
      "source_no_pretreatment_flower_cover_assessed":True,
      "source_no_plant_compatible_pollen_or_seed_measured":True,
      "source_dataset_has_unmatched_species_arm_selection":True,
      "results":result,
      "wet_macro_recording_bias_mechanism_identified":False
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    j=analyse(fetch(),None)
    out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(j,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    report={}
    for t,x in j["results"].items():
        report[t]={
           "plant_source_n":x["original_source_scored_matched_rows"],
           "wet":x["original_wet_plant_names"],
           "dry":x["original_drier_plant_names"],
           "same_meadow_same_event":x["within_meadow_event_matched_wet_drier_flower"],
           "single_arm_selection":x["original_single_arm_plant_recording_selection"],
           "dry_only":x["dry_only_F_le7_flower_association"],
           "by_meadow":x["per_physical_meadow_within_event_source_species_support"],
        }
    print("GIAVI2021_ORIGINAL_SAME_MEADOW_EVENT_WET_DRY_FLORAL_NULL_AUDIT "+
          json.dumps(report,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
