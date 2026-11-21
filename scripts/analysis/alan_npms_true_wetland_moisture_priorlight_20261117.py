#!/usr/bin/env python3
"""Prespecified actual NPMS-site habitat × plant MOIST × prior PCNL.

Use only original 2015–25 NPMS habitat labels (same survey/source period),
exact official European EIV and the previously verified 61,555 record positive
cover source. Infer no species absence, plant fitness or direct lighting effect.
The habitat label is stable within a physical NPMS plot across years or the
plot is excluded. Inland wetland categories are named before opening outcome.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from scripts.acquisition import alan_npms_original_habitat_source_gate_20261117 as h
from scripts.analysis import alan_npms_independent_moisture_priorlight_20261117 as core

SCHEMA="alan_npms_measured_habitat_wetness_light_candidate_v1"
INLAND_WET=("Freshwater","Marsh and fen","Bog and wet heath")
EXCLUDE=("Not in scheme","Coast")
MIN_WET_PLOTS=30
MIN_REF_PLOTS=100
MIN_WET_MIXED=80
MIN_REF_MIXED=200
MIN_WET_SPECIES=15
MIN_DRIER_SPECIES=15

def study(npms,cover,receipt,eiv,return_source_panels=False):
    source,zfull=h.analyze(npms,cover,receipt,return_matched=True)
    data,source_cover=core.load(cover,receipt,eiv)
    # Source original positive-cover row identity is one physical plot ×
    # source survey × season × species × year; do not join nearest locality.
    k=["location_id","survey_method_id","season","species","year"]
    pairs=zfull[k+["broad"]].copy()
    assert not pairs.duplicated(k).any(),"habitat source duplicates botanical plot rows"
    # Source landcover membership is defined with all source-positive records
    # and kept fixed over all years, NOT selected using plant cover outcomes.
    site_labels=zfull.groupby("location_id").broad.agg(
        lambda s:sorted(set(v for v in s if v)))
    stable=site_labels.loc[site_labels.apply(len)==1].map(lambda x:x[0])
    allowed=stable[~stable.isin(EXCLUDE)]
    joins=data.merge(pairs,on=k,how="left",validate="one_to_one",indicator=True)
    assert len(joins)==source_cover["moisture_source_exact_rows"]
    joins["site_class"]=joins.location_id.map(allowed)
    # A physical plot may have one stable named habitat among documented years
    # but some other recorded survey years have UNKNOWN habitat. Never fill
    # those rows from adjacent years: keep only exact same-year habitat labels.
    eligible=joins.loc[joins.site_class.notna() & joins.broad.eq(joins.site_class)].copy()
    assert eligible.broad.eq(eligible.site_class).all(),"original survey habitat vs stable site discrepancy"
    eligible["inland_wetland"]=eligible.site_class.isin(INLAND_WET)
    eligible["light_moist_wetland"]=eligible.light_moist*eligible.inland_wetland.astype(int)
    # Avoid habitat × light source confounding by physically matched same
    # plot/survey/season/year wet + drier plant observations.
    siteyears=eligible.groupby("plotyear_method_season").agg(
        wet=("wet7","max"),low=("wet7","min"),
        plot=("location_id","first"),inland=("inland_wetland","first"),
        unique_plants=("species","nunique"))
    mixed=siteyears[(siteyears.wet==1)&(siteyears.low==0)&(siteyears.unique_plants>=2)]
    common=eligible.loc[eligible.plotyear_method_season.isin(mixed.index)].copy()
    def support(d):
        return {"positive_cover_rows":len(d),
                "physical_plots":int(d.location_id.nunique()),
                "source_plant_species":int(d.species.nunique()),
                "plant_moist_gt7_species":int(d.loc[d.wet7==1,"species"].nunique()),
                "plant_moist_le7_species":int(d.loc[d.wet7==0,"species"].nunique()),
                "plotyear_methodseason_with_both_plant_moist_groups":int(
                    d.plotyear_method_season.nunique())}
    sample={
      "schema":SCHEMA,
      "source_npms_CRC_and_positive_cover_joins_verified":True,
      "source_years":[2015,2025],
      "inland_wetland_original_categories":list(INLAND_WET),
      "excluded_original_categories":list(EXCLUDE),
      "source_original_NPMS_habitat_scored_plots":source["source_cover_unique_broad_habitat_plot_count"],
      "true_original_positive_cover_before_moisture":source["original_cover_rows"],
      "source_original_cover_scored":source_cover,
      "stable_same_original_habitat_sites":len(stable),
      "stable_accepted_category_sites":len(allowed),
      "eligible_all_known_species":support(eligible),
      "wetland_known_species":support(eligible[eligible.inland_wetland]),
      "other_habitat_known_species":support(eligible[~eligible.inland_wetland]),
      "wetland_both_guilds":support(common[common.inland_wetland]),
      "other_habitat_both_guilds":support(common[~common.inland_wetland]),
      "original_illuminated_treatment_randomized":False,
      "species_not_reported_was_never_assigned_absence":True,
      "estimated_actual_reproductive_or_visitation_response":False
    }
    wet=sample["wetland_both_guilds"];dry=sample["other_habitat_both_guilds"]
    sample["source_gate"]={
        "min_wetland_plots":MIN_WET_PLOTS,
        "min_terrestrial_plots":MIN_REF_PLOTS,
        "min_wetland_mixed_plotyear":MIN_WET_MIXED,
        "min_terrestrial_mixed_plotyear":MIN_REF_MIXED,
        "min_each_moisture_guild_species":MIN_WET_SPECIES
    }
    supported=(wet["physical_plots"]>=MIN_WET_PLOTS and dry["physical_plots"]>=MIN_REF_PLOTS and
         wet["plotyear_methodseason_with_both_plant_moist_groups"]>=MIN_WET_MIXED and
         dry["plotyear_methodseason_with_both_plant_moist_groups"]>=MIN_REF_MIXED and
         all(s["plant_moist_gt7_species"]>=MIN_WET_SPECIES and
             s["plant_moist_le7_species"]>=MIN_DRIER_SPECIES for s in (wet,dry)))
    sample["source_gate"]["pass"]=bool(supported)
    if not supported:
        sample["status"]="source_insufficient_wetland_or_reference_common_support_stop"
        return sample
    if return_source_panels:
        # Internal, in-memory reuse only; source records must never be
        # committed or uploaded. The original survey IDs remain protected.
        return sample,common,eligible
    # The common plot-year FE absorbs both site-year light level and any
    # common hydro-climate/observer/season signal. The residual light-MOIST
    # slope is descriptive and not a causal biological ALAN effect.
    F=["panel_key","species_year","plotyear_method_season"]
    sample["model_source_rows"]=len(common)
    sample["model_wetland"]=core.fit(
        common[common.inland_wetland],["light_moist","doy30"],F,two_way=True)
    sample["model_other_habitats"]=core.fit(
        common[~common.inland_wetland],["light_moist","doy30"],F,two_way=True)
    sample["model_three_way_delta"]=core.fit(
        common,["light_moist","light_moist_wetland","doy30"],F,two_way=True)
    # POSTHOC source-strength and estimator sensitivity, explicitly report
    # all signs and STOPs. Thresholds .1/.3 are inherited from the older
    # independent NPMS analysis, not chosen on wetland coefficient signs.
    sample["posthoc_robustness"]={}
    sample["posthoc_robustness"]["record_equal"]=core.fit(
        common,["light_moist","light_moist_wetland","doy30"],F,"record_equal",two_way=True)
    panel=eligible.groupby("panel_key").prior_light.agg(["min","max"])
    for cutoff in (.1,.3):
        keep_ids=panel.index[(panel["max"]-panel["min"])>=cutoff]
        pool=eligible.loc[eligible.panel_key.isin(keep_ids)]
        cells=pool.groupby("plotyear_method_season").agg(
            highest=("wet7","max"),lowest=("wet7","min"),n=("species","nunique"))
        both=cells.index[(cells.highest==1)&(cells.lowest==0)&(cells.n>=2)]
        pool=pool.loc[pool.plotyear_method_season.isin(both)].copy()
        w=support(pool[pool.inland_wetland]);o=support(pool[~pool.inland_wetland])
        gate=(w["physical_plots"]>=30 and o["physical_plots"]>=100 and
              w["plotyear_methodseason_with_both_plant_moist_groups"]>=65 and
              o["plotyear_methodseason_with_both_plant_moist_groups"]>=150 and
              min(w["plant_moist_gt7_species"],w["plant_moist_le7_species"],
                  o["plant_moist_gt7_species"],o["plant_moist_le7_species"])>=12)
        item={"source_threshold":cutoff,"inland_wetland":w,"other_habitat":o,
              "eligibility_passed":bool(gate),"posthoc":True}
        if gate:item["fit"]=core.fit(
            pool,["light_moist","light_moist_wetland","doy30"],F)
        else:item["fit"]={"status":"not_enough_original_source_common_support_STOP"}
        sample["posthoc_robustness"][str(cutoff)]=item
    # Category-leverage: these are mutually overlapping sensitivity fits,
    # not independent experiments. Any failure to support exclusion is STOP,
    # never rescue by lowering sample criteria on the effect sign.
    sample["posthoc_robustness"]["leave_one_inland_class_out"]={}
    for excluded in INLAND_WET:
        pool=common.loc[~(
            common.inland_wetland & common.site_class.eq(excluded))].copy()
        w=support(pool[pool.inland_wetland])
        other=support(pool[~pool.inland_wetland])
        gate=(w["physical_plots"]>=25 and
              w["plotyear_methodseason_with_both_plant_moist_groups"]>=60 and
              w["positive_cover_rows"]>=500 and
              min(w["plant_moist_gt7_species"],w["plant_moist_le7_species"])>=12 and
              other["physical_plots"]>=100)
        item={"remaining_wet":w,"source_gate":bool(gate),"posthoc":True}
        if gate:item["fit"]=core.fit(
            pool,["light_moist","light_moist_wetland","doy30"],F)
        else:item["fit"]={"status":"cannot_leave_class_out_with_adequate_source_STOP"}
        sample["posthoc_robustness"]["leave_one_inland_class_out"][excluded]=item
    sample["status"]="source_supported_exploratory_positive_cover_original_habitat_moderator"
    return sample

def main():
    a=argparse.ArgumentParser()
    for k in ("npms","cover","receipt","eiv","out"):
        a.add_argument("--"+k,required=True)
    q=a.parse_args()
    s=study(q.npms,q.cover,q.receipt,q.eiv)
    p=Path(q.out);p.parent.mkdir(exist_ok=True,parents=True)
    p.write_text(json.dumps(s,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_NPMS_HABITAT_WETNESS_LIGHT_GATES",json.dumps({
        "status":s["status"],"support_wet":s["wetland_both_guilds"],
        "support_other":s["other_habitat_both_guilds"],
        "threeway":s.get("model_three_way_delta",{}).get("coefficients",{}),
        "slope_wet":s.get("model_wetland",{}).get("coefficients",{}),
        "slope_other":s.get("model_other_habitats",{}).get("coefficients",{})
    },sort_keys=True))
if __name__=="__main__":main()
