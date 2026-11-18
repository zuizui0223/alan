#!/usr/bin/env python3
"""2026 Lu original 45 plant-label network vs 11 focal-species table, SOURCE-ONLY.

Retains precisely the 1,605 source rows with complete treatment/diel IDs
from the original 1,606-row 2026 Data S1 "plant–pollinator network" sheet.
The one source row with missing ALAN status and diel must be excluded, NEVER
forward-filled from adjacent rows. All valid records preserve original
plant×spatial block×survey-round×light×diel; missing plant ≠ zero visits.

Source original workbook sha256:
ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f

Do not pool with author "visitation" 400-row source or selected focal 1038-row
species table, whose raw Diptera sums can differ, although the five spatial
light treatment blocks are the same. This audit estimates visitation
redistribution, NOT stigma pollen receipt, selfing, reproduction or wetland
plant abundance. Its many plants are NOT independent ALAN experiments.
"""
from __future__ import annotations
import argparse
import collections
import json
import math
import re
from pathlib import Path

from scripts.acquisition.alan_lu2026_original_workbook_schema import source_open,sheet_data,ORIGINAL_WORKBOOK_SHA
from scripts.analysis.alan_lu2026_original_fiveblock_guild_replacement import (
    rows,num,normalize_time,normalize_light,need
)
from scripts.analysis.alan_lu2026_original_sameplant_guild_replacement import (
    source_results,read_author_taxon_abbreviations
)

SCHEMA="alan_lu2026_original_45plant_network_same_flower_20261118_v1"
NETWORK_SHEET="plant–pollinator network"
FIELDS=("all","Lepidoptera","Diptera","Hymenoptera","Coleoptera","Others")
SOURCE_ROWS=1606
COMPLETE_RECORDS=1605

def parse_source(raw):
    src=rows(raw)
    need(len(src)==SOURCE_ROWS,"original 1606-row full plant network changed")
    required={"plant","time","sampling_time","treatment","block","repeat",*FIELDS}
    need(required.issubset(src[0]),"original full network source columns changed")
    out=[]
    excluded=collections.Counter()
    sum_guild_discrepancies=0
    for row in src:
        phase=normalize_time(row["sampling_time"])
        arm=normalize_light(row["treatment"])
        plant=str(row["plant"]).strip()
        block=str(row["block"]).strip()
        period=str(row["time"]).strip()
        repeat=str(row["repeat"]).strip()
        # Author original contains ONE row with no source lighting/diel
        # coding. No information is carried down from previous row.
        if phase.startswith("unknown:") or arm.startswith("unknown:"):
            excluded["source_row_unknown_phase_or_treatment"]+=1
            continue
        need(phase in ("day","night") and arm in ("lit","dark"),
             "unsupported original diel/treatment value")
        need(plant and block and period and repeat,
             "valid original source missing plant/spatial block/time")
        counts={k:num(row[k]) for k in FIELDS}
        need(all(x is not None and x>=0 and x.is_integer() for x in counts.values()),
             "original valid source guild count is missing/noninteger")
        if counts["all"]!=sum(counts[k] for k in FIELDS if k!="all"):
            sum_guild_discrepancies+=1
        out.append({"plant":plant,"phase":phase,"arm":arm,"block":block,
                    "round":period,"repeat":repeat,**counts})
    need(len(out)==COMPLETE_RECORDS and excluded["source_row_unknown_phase_or_treatment"]==1,
         "original network missing source keys changed; do not impute")
    need(sum_guild_discrepancies==0,"original guild subtotals do not equal original row total")
    blocks={x["block"] for x in out}
    need(len(blocks)==5,"source five independent lighting blocks changed")
    return out,{"source_original_worksheet_raw_rows":len(src),
                 "source_valid_plant_visitation_rows":len(out),
                 "original_rows_with_unreported_light_or_diel":excluded[
                     "source_row_unknown_phase_or_treatment"],
                 "original_source_guild_subtotal_row_disagreements":sum_guild_discrepancies,
                 "original_distinct_plant_labels":len({r["plant"] for r in out}),
                 "original_plant_labels_by_diel":{
                     x:len({r["plant"] for r in out if r["phase"]==x})
                     for x in ("night","day")},
                 "source_plants_are_not_independent_lighting_installations":True,
                 "original_unlabelled_plant_not_assumed_0_visits":True}

def compare_views(full,visit,focal):
    def totals(rows0,field_map):
        x=collections.Counter()
        for v in rows0:
            phase=normalize_time(v["sampling_time"])
            if phase not in ("night","day"):continue
            for guild,field in field_map.items():
                a=num(v[field])
                if a is not None:x[(phase,guild)]+=a
        return x
    source_total=totals(rows(visit),{
        "all":"abundance","Diptera":"D_abundance",
        "Lepidoptera":"L_abundance","Hymenoptera":"H_abundance",
        "Coleoptera":"C_abundance","Others":"O_abundance"})
    focal_total=totals(rows(focal),{
        "all":"all","Diptera":"Diptera","Lepidoptera":"Lepidoptera",
        "Hymenoptera":"Hymenoptera","Coleoptera":"Coleoptera","Others":"Others"})
    full_total=collections.Counter()
    for r in full:
        for g in FIELDS:full_total[(r["phase"],g)]+=r[g]
    result={}
    for stage in ("night","day"):
        result[stage]={
            g:{"full_1605_row_network":full_total[(stage,g)],
               "original_400_row_visitation_summary":source_total[(stage,g)],
               "original_1038_focal_species_visitation":focal_total[(stage,g)],
               "full_network_minus_summary":full_total[(stage,g)]-source_total[(stage,g)],
               "focal_minus_summary":focal_total[(stage,g)]-source_total[(stage,g)]}
            for g in FIELDS}
    return {
        "source_three_views_exact_same_visitation_totals":all(
           full_total[(s,g)]==source_total[(s,g)]==focal_total[(s,g)]
           for s in ("night","day") for g in FIELDS),
        "source_plant_and_guild_count_per_view":result,
        "do_not_sum_across_three_views":True
    }

def execute():
    workbook,rcode,src=source_open()
    need(src["source_workbook_sha256"]==ORIGINAL_WORKBOOK_SHA,
         "published Lu workbook changed from source SHA")
    sheets=sheet_data(workbook)
    full,support=parse_source(sheets[NETWORK_SHEET])
    author_map=read_author_taxon_abbreviations(
        rcode,{row["plant"] for row in full})
    result=source_results(full,sorted({x["block"] for x in full}),author_map)
    original=compare_views(full,sheets["visitation"],sheets["species visitation"])
    # Same-light treatment in each of five independent blocks, not
    # interpreting within-source subrows as independent interventions.
    return {"schema":SCHEMA,"source_original_DOI":"10.1002/ece3.74368",
            "original_data_sha256":ORIGINAL_WORKBOOK_SHA,
            "source_spatial_pairs":5,"source_original_measured_table":support,
            "original_author_R_explicit_name_mapping_only":author_map,
            "same_plant_study_within_network_view":result,
            "three_views_denominator_audit":original,
            "moth_loss_fly_gain_is_visitation_not_functional_pollen":True,
            "source_original_species_seed_set_not_measured":True,
            "undocumented_species_or_missing_light_arm_not_assigned_zero":True,
            "meta_independent_species_replication_not_established":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    result=execute()
    dest=Path(args.out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("LU2026_FULL_ORIGINAL_NETWORK_SAME_FLOWER_GUILD_REDISTRIBUTION",json.dumps({
        "original_source":result["source_original_measured_table"],
        "source_denominators_identical":
            result["three_views_denominator_audit"]["source_three_views_exact_same_visitation_totals"],
        "coincident_plant_coverage":{k:{
          "source_plant_labels":v["original_distinct_matched_plants"],
          "night_plant_block_pairs":v["source_plants_same_block_treatment_pairs_by_diel"]["night"],
          "night_moth_down_fly_up":v["night_same_plant_moth_down_fly_up"],
          "strongest_original_source_plants":v["night_largest_documented_species_changes_descriptive"][:12]}
          for k,v in result["same_plant_study_within_network_view"].items()},
        "visitation_totals_comparison_night":
            result["three_views_denominator_audit"]["source_plant_and_guild_count_per_view"]["night"]
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
