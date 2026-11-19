#!/usr/bin/env python3
"""Six primary programme A evidence map, NOT a derived global wetland meta.

Phenomenological reversal across 2017-2025 Westhavelland papers is
NOT a single same-taxon effect reversal: one lit and one dark site, studies
have different years and emergence/flying endpoints.

Water-land moderators require >1 independent source programme and
matched response/process. This map deliberately refuses to calculate
a meta or the forest vs wetland moderator from these six case studies.
"""
import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path

SOURCE=Path("config/alan_stageA_water_land_process_specific_primary_gate_20261119.json")

def audit(data):
    if data.get("schema")!="alan_stageA_primary_original_water_land_outcome_and_programme_partition_v1":
        raise ValueError("Original stageA process source schema drift")
    progs=data.get("programmes",[])
    pids=[r.get("programme_id") for r in progs]
    if len(pids)!=6 or len(set(pids))!=6:
        raise ValueError("Six source-targeted original experimental programme IDs not unique")
    expected={
      "Parkinson2020_20freshwater_mesocosms",
      "Perkin2014_4paired_headwater_streams",
      "WesthavellandDrainageDitch2012",
      "NavarroBarranco2015_GreatBarrierReef_lighttraps",
      "Deitsch2023_hardwood_4pairs",
      "Farnworth2026_NZ_forest_predators",
    }
    if set(pids)!=expected:
        raise ValueError("Original source mechanistic programme identities changed")
    west=next(r for r in progs if r["programme_id"]=="WesthavellandDrainageDitch2012")
    if (west["original_n_treated_sites"]!=1 or west["original_n_dark_sites"]!=1
        or set(west["same_programme_source_papers"])!={
            "10.3389/fenvs.2017.00061","10.1007/s00027-025-01161-7"}):
        raise ValueError("Westhavelland is one lit and one dark site, not two experiments")
    speciesstage={x["stage"]:x for x in west["source_summary"]}
    emergence=speciesstage["Diptera_emergence_2014"]
    ratio=emergence["light_dark_mean_lit"]/emergence["light_dark_mean_dark"]
    if abs(ratio-3.72/7.84)>1e-10 or not ratio>0:
        raise ValueError("Original aquatic Diptera 2025 lit/dark mean drift")
    perkin=next(r for r in progs if r["programme_id"]=="Perkin2014_4paired_headwater_streams")
    if (perkin["original_n_independent_streams"]!=4 or
        perkin["original_n_light_reach_per_stream"]!=1 or
        perkin["original_n_dark_reach_per_stream"]!=1):
        raise ValueError("Original four paired streams light source is nested within streams")
    parkinson=next(r for r in progs if r["programme_id"]=="Parkinson2020_20freshwater_mesocosms")
    if (parkinson["original_n_light_installation"]!=10 or
        parkinson["original_n_dark_installation"]!=10):
        raise ValueError("Original freshwater mesocosm 10+10 design invalid")
    marine=next(r for r in progs if r["programme_id"].startswith("NavarroBarranco"))
    if marine["site_type"]!="marine_shallow_water":
        raise ValueError("Do not pool marine amphipod light-trap mobility with freshwater aquatic emergence")
    habitats=Counter(r["site_type"] for r in progs)
    direct_process=defaultdict(set)
    for r in progs:
        for outcome in r["source_summary"]:
            direct_process[outcome["stage"]].add(r["programme_id"])
    source_doi=set()
    for r in progs:
        source_doi.update(r.get("same_programme_source_papers") or [r["doi"]])
    if len(source_doi)!=7:
        raise ValueError("Six source programmes expected seven scientific publications")
    gates=data["source_inference_gates"]
    if not all(gates.values()):
        raise ValueError("Original phase A response independence inference gates changed")
    return {
      "schema":"alan_20261119_water_land_six_source_programmes_distinct_process_no_meta_v1",
      "targeted_original_programme_k":6,
      "targeted_publication_n":len(source_doi),
      "original_programme_habitat_k_not_a_meta_denominator":dict(sorted(habitats.items())),
      "original_programmes_per_phenomenological_measurement_NOT_comparable_effect":{
          x:len(uu) for x,uu in sorted(direct_process.items())
      },
      "Westhavelland_moisture_context":"One lit site and one dark site; 2017 and 2025 source papers same physical infrastructure",
      "Westhavelland_2014_Diptera_emergence_CPUE_lit_dark_descriptive_mean_ratio":ratio,
      "Westhavelland_2014_Diptera_emergence_descriptive_log_mean_ratio_NO_INDEPENDENT_SITE_SE":math.log(ratio),
      "Westhavelland_2014_flying_Diptera_catch_approximately_light_dark_ratio":8,
      "Westhavelland_prior_aquatic_origin_2012_2013_riparian_emergence":"increase reported but not same Diptera-only/year measurement as 2025",
      "Perkin2014_4_stream_pairs_night_drift_rough_lit_dark_ratio":.5,
      "Perkin2014_emergence_and_trout_growth":"not resolved in published source, NOT zero",
      "Parkinson2020_20_mesocosms_web_spider_ratio":2.01,
      "Parkinson2020_20_mesocosms_aquatic_origin_flying_insect_ratio":6.63,
      "marine_amphipod_lighttrap_attraction_not_freshwater_emergence":True,
      "independent_A_wetland_by_terrestrial_meta_moderator_fit":False,
      "study_targeted_not_PRISMA_complete":True,
      "no_sample_size_from_nested_streetlamps_or_traps":True,
      "no_meta_pooled_mean_CI_or_pvalue":True,
      "no_habitat_sign_claim_without_matched_endpoint_and_independent_site_replicates":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SOURCE)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    d=audit(json.loads(args.source.read_text()))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_STAGE_A_SIX_PROGRAMME_HABITATS",d["original_programme_habitat_k_not_a_meta_denominator"])
    print("ALAN_STAGE_A_WESTHAVELLAND_TRUE_RIPARIAN_DIPTERA_EMERGENCE_RATIO",
          d["Westhavelland_2014_Diptera_emergence_CPUE_lit_dark_descriptive_mean_ratio"])
    print("ALAN_STAGE_A_NO_PSEUDOREPLICATION_META_GATE",d["independent_A_wetland_by_terrestrial_meta_moderator_fit"])


if __name__=="__main__":
    main()
