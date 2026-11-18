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
    source_results,read_author_taxon_abbreviations,match_plants
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

def theoretical_same_flower_fly_effectiveness_threshold(full):
    """Hypothetical q_fly/q_moth needed to replace moth-visit loss.

    Under an explicitly counterfactual linear same-flower pollen model
    ΔP = q_moth*Δmoth + q_fly*Δfly, assuming q_moth>0, stable per-visit
    effectiveness across light conditions, no other guild changes,
    matched survey effort, and no autonomous selfing contribution change.
    This is NOT an observed pollen deposition or seed-set quantity.
    """
    paired=match_plants(full,[])["same_plant_block_round_repeat"]
    byplant=collections.defaultdict(lambda:{
        "blocks":set(),"moth_delta":0.,"fly_delta":0.,
        "all_delta":0.,"source_matched_periods":0})
    for r in paired:
        if r["phase"]!="night":continue
        x=byplant[r["plant"]]
        x["blocks"].add(r["block"])
        x["moth_delta"]+=r["moth_delta"]
        x["fly_delta"]+=r["fly_delta"]
        x["all_delta"]+=r["total_delta"]
        x["source_matched_periods"]+=1
    eligible={}
    for plant,x in byplant.items():
        nb=len(x["blocks"])
        if nb<2:continue
        ml=x["moth_delta"]
        fg=x["fly_delta"]
        # If flies do not increase, no positive finite fly effectiveness
        # can directly compensate a lost moth under this two-guild scenario.
        threshold=(-ml/fg) if ml<0 and fg>0 else None
        eligible[plant]={
            "observed_same_lit_dark_spatial_blocks":nb,
            "observed_paired_survey_periods":x["source_matched_periods"],
            "source_moth_visits_lit_minus_dark":ml,
            "source_fly_visits_lit_minus_dark":fg,
            "source_total_visits_lit_minus_dark":x["all_delta"],
            "counterfactual_min_q_fly_per_visit_relative_to_q_moth":threshold,
            "threshold_defined_only_if_moth_lost_and_fly_gained":ml<0 and fg>0,
            "no_direct_pollen_deposition_or_seed_set_observed":True}
    defined=[v for v in eligible.values()
             if v["counterfactual_min_q_fly_per_visit_relative_to_q_moth"] is not None]
    return {
       "original_species_with_2plus_independent_light_blocks":len(eligible),
       "species_with_moth_loss_and_fly_gain_under_same_survey_round":len(defined),
       "source_original_plant_specific_counterfactual_thresholds":eligible,
       "equation":"ΔP = q_moth*ΔV_moth + q_fly*ΔV_fly; if ΔV_moth<0 and ΔV_fly>0 then q_fly/q_moth>=-ΔV_moth/ΔV_fly",
       "assumptions":"Linear no-saturation effective pollen deposition; q constant across light/dark, other guilds unchanged, autonomous selfing unchanged, matched exposure; NONE directly verified",
       "not_measurement_of_fly_or_moth_pollen_efficiency":True,
       "denominator_is_original_same_plant_same_block_same_round_visits":True,
       "no_species_level_hypothesis_testing_or_plant_fitness_prediction":True
    }

def joint_day_night_same_plant_routes(full):
    """Can daytime or nighttime Diptera visit the SAME plant losing night moths?"""
    paired=match_plants(full,[])
    all_by_pair={}
    for row in paired["same_plant_block_any_round"]:
        key=(row["plant"],row["block"])
        d=all_by_pair.setdefault(key,{})
        # Requiring the SAME plant at SAME block with original lit/dark
        # recordings for each diel session separately.
        if row["phase"] in d:raise ValueError("original same plant×block×phase duplicated")
        d[row["phase"]]=row
    checked=[(plant,block,v["night"],v["day"])
             for (plant,block),v in all_by_pair.items()
             if "night" in v and "day" in v]
    stat=collections.Counter()
    plants=collections.defaultdict(set)
    spatial_blocks=collections.defaultdict(set)
    original_night_moth_loss=0.0
    original_day_fly_gain=0.0
    original_night_fly_gain=0.0
    documented=[]
    for plant,block,night,day in checked:
        moth_down=night["moth_delta"]<0
        noct_fly_up=night["fly_delta"]>0
        diurnal_fly_up=day["fly_delta"]>0
        categories={
           "night_moth_down":moth_down,
           "night_moth_down_and_night_fly_up":moth_down and noct_fly_up,
           "night_moth_down_and_day_fly_up":moth_down and diurnal_fly_up,
           "night_moth_down_and_either_fly_up":
               moth_down and (noct_fly_up or diurnal_fly_up),
           "night_moth_down_and_both_fly_up":
               moth_down and noct_fly_up and diurnal_fly_up,
           "night_moth_down_without_fly_increase":
               moth_down and not noct_fly_up and not diurnal_fly_up
        }
        for cat,yes in categories.items():
            if yes:
                stat[cat]+=1
                plants[cat].add(plant)
                spatial_blocks[cat].add(block)
        if moth_down:
            original_night_moth_loss-=night["moth_delta"]
            if noct_fly_up:original_night_fly_gain+=night["fly_delta"]
            if diurnal_fly_up:original_day_fly_gain+=day["fly_delta"]
        documented.append((plant,block,night["moth_delta"],night["fly_delta"],
                           day["fly_delta"],noct_fly_up,diurnal_fly_up))
    return {
      "same_plant_and_spatial_block_has_both_arms_for_day_and_night":len(checked),
      "distinct_original_plant_labels":len({k for k,b,_,_ in checked}),
      "distinct_independent_experimental_light_blocks":len({b for k,b,_,_ in checked}),
      "patterns":{cat:{
            "source_plant_block_pairs":stat[cat],
            "distinct_plant_labels":len(plants[cat]),
            "independent_light_blocks_with_pattern":len(spatial_blocks[cat])
          } for cat in (
           "night_moth_down",
           "night_moth_down_and_night_fly_up",
           "night_moth_down_and_day_fly_up",
           "night_moth_down_and_either_fly_up",
           "night_moth_down_and_both_fly_up",
           "night_moth_down_without_fly_increase")},
      "total_night_moth_visit_decline_among_night_moth_declining_same_plant_pairs":
          original_night_moth_loss,
      "positive_night_fly_visit_gain_on_those_same_plant_pairs":
          original_night_fly_gain,
      "positive_day_fly_visit_gain_on_those_same_plant_pairs":
          original_day_fly_gain,
      "different_diel_visit_count_denominators_not_interchangeable":True,
      "same_flower_temporal_visitor_change_not_a_pollen_or_seed_rescue":True,
      "five_light_blocks_not_plant_block_independent_studies":True
    }

def within_fullnetwork_five_block_guilds(full):
    # Counts within the 1,605 valid *same full network worksheet* rows.
    # There are only FIVE independent original spatial lit/dark pairs.
    from scripts.analysis.alan_lu2026_original_fiveblock_guild_replacement import exact_five_signflip
    blocks=sorted({x["block"] for x in full})
    need(len(blocks)==5,"full author network no longer has five blocks")
    by=collections.Counter()
    for x in full:
        for guild in ("all","Lepidoptera","Diptera"):
            by[(x["phase"],x["arm"],x["block"],guild)]+=x[guild]
    out={}
    for diel in ("night","day"):
        guilds={}
        for guild in ("all","Lepidoptera","Diptera"):
            lit=[by[(diel,"lit",b,guild)] for b in blocks]
            dark=[by[(diel,"dark",b,guild)] for b in blocks]
            delta=[math.log((a+0.5)/(c+0.5)) for a,c in zip(lit,dark)]
            x=exact_five_signflip(delta)
            guilds[guild]={
              "sum_lit_visits":sum(lit),
              "sum_dark_visits":sum(dark),
              "lit_dark_total_visits_ratio":
                 sum(lit)/sum(dark) if sum(dark)>0 else None,
              "five_pair_geometric_lit_dark_visits_ratio":
                 math.exp(x["paired_5_block_mean"]),
              "descriptive_two_sided_32_signflip_p":
                 x["descriptive_exact_32_paired_label_signflips_two_sided_p"],
              "block_count":len(blocks),
              "original_experimental_light_units_not_plant_rows":10
            }
        moth=[math.log((by[(diel,"lit",b,"Lepidoptera")]+.5)/
                       (by[(diel,"dark",b,"Lepidoptera")]+.5)) for b in blocks]
        flies=[math.log((by[(diel,"lit",b,"Diptera")]+.5)/
                        (by[(diel,"dark",b,"Diptera")]+.5)) for b in blocks]
        joint=exact_five_signflip([a-b for a,b in zip(moth,flies)])
        out[diel]={
           "original_author_full_network_source_guilds":guilds,
           "five_block_relative_moth_to_fly_lit_vs_dark_geometric_ratio":
               math.exp(joint["paired_5_block_mean"]),
           "five_block_moth_fly_relative_signflip_p":
               joint["descriptive_exact_32_paired_label_signflips_two_sided_p"],
           "relative_moth_to_fly_direction_negative_blocks":
               joint["independent_block_negative_effects"],
           "visitation_not_functional_pollen_substitution":True
        }
    return out

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
            "source_temporal_day_night_same_flower_substitution":
                joint_day_night_same_plant_routes(full),
            "hypothetical_pollen_visit_effectiveness_thresholds":
                theoretical_same_flower_fly_effectiveness_threshold(full),
            "original_full_network_only_five_block_guilds":
                within_fullnetwork_five_block_guilds(full),
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
        "day_night_same_flower_routes":
            result["source_temporal_day_night_same_flower_substitution"],
        "potential_fly_effectiveness_thresholds":
            result["hypothetical_pollen_visit_effectiveness_thresholds"],
        "original_full_network_night_guild_shift":
            result["original_full_network_only_five_block_guilds"]["night"],
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
