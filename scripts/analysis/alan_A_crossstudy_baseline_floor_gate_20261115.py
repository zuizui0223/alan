#!/usr/bin/env python3
"""Cross-study ALAN reproduction hurdle hypothesis audit (two preprints, no meta).

Li: re-fetch author-signed 79 plants, 11 matched three-arm L/P/C strata.
Grognuz: local 28 source-checked supplementary blocks (not raw 448 pots).
Dark-control zeros produce nonnegative lit-dark differences by construction.
Groups are not independent lighting treatments and cannot be pooled.
"""
from __future__ import annotations
import json,math
from pathlib import Path
from scripts.analysis import alan_Li2026_original_triple_arm_20261115 as li
from scripts.analysis import alan_Grognuz2026_block_floor_signs_20261115 as grog
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"results/alan_A_crossstudy_baseline_zero_hurdle_gate_20261115.json"
SOURCE_REGISTER=ROOT/"config/alan_A_reproductive_screen_20261115.json"

def verify(liout,grout,registry):
    ids={p["id"] for p in registry["entries"] if p["light_design"].startswith("manipulated")
         and p["a_primary_stage"]!="exclude_A_no_reproduction"}
    if ids!={"knop2017","macgregor2019","boom2020","giavi2020","wilson2021",
            "li2026","lu2026jpe","xiao2025","grognuz2026pre"}:
        raise AssertionError("Frozen nine experimental papers drifted")
    a=liout["all_three_arm_matched"]
    t=grout["Trifolium_pratense"];c=grout["Centaurea_jacea"]
    if liout["source"]["original_plants"]!=79 or liout["source"]["original_seeds_observed"]!=78:
        raise AssertionError("Li original source drift")
    if liout["source"]["three_arm_strata"]!=11 or a["n_dark_zero"]!=5:
        raise AssertionError("Li true 3-arm matching drift")
    if a["seeds_LminusP"]["strata_signs"]!={"positive":8,"negative":3,"zero":0}:
        raise AssertionError("Li positive/negative signed treatment comparisons drift")
    if a["seeds_LminusP"]["field_signs"]!={"positive":2,"negative":1,"zero":0}:
        raise AssertionError("Li three true field signs drift")
    if grout["total_blocks"]!=28 or t["n_blocks"]!=20 or c["n_blocks"]!=8:
        raise AssertionError("Grognuz source 28 choice blocks drift")
    if t["all_block_signs"]!={"lit_positive":8,"lit_negative":12,"zero":0}:
        raise AssertionError("Trifolium source block sign drift")
    if c["all_block_signs"]!={"lit_positive":3,"lit_negative":5,"zero":0}:
        raise AssertionError("Centaurea source block sign drift")
    if len(t["dark_zero_blocks"])!=1 or len(c["dark_zero_blocks"])!=1:
        raise AssertionError("Grognuz zero-dark groups changed")
    return True

def result(liout,grout,reg):
    verify(liout,grout,reg)
    la=liout["all_three_arm_matched"]
    t=grout["Trifolium_pratense"];c=grout["Centaurea_jacea"]
    summary={
      "schema":"ALAN_A_20261115_two_experiment_source_signed_baseline_floor_gate_v1",
      "study_publications":2,
      "distinct_designs":2,
      "confirmed_manipulated_reproduction_census_unchanged":9,
      "Li2026":{
       "source_git_blob":liout["author_git_blob"],
       "original_plant_rows":liout["source"]["original_plants"],
       "original_seed_outcomes_known":liout["source"]["original_seeds_observed"],
       "true_independent_fields":3,
       "field_taxon_round_triplets":la["n"],
       "control_seed_zero_triplets":la["n_dark_zero"],
       "seeds_per_initial_bud_ratios":{
         "linear_over_dark":la["seeds_LoverC_pooled_rate_ratio"],
         "point_over_dark":la["seeds_PoverC_pooled_rate_ratio"],
         "linear_over_point":la["seeds_LoverP_pooled_rate_ratio"]},
       "stratum_signs_L_over_P":la["seeds_LminusP"]["strata_signs"],
       "true_field_signs_L_over_P":la["seeds_LminusP"]["field_signs"],
       "minimum_two_plants_each_treatment":{
           "source_strata":liout["at_least_two_plants_per_each_L_P_C_arm"]["n"],
           "L_over_P_seed_stratum_signs":
             liout["at_least_two_plants_per_each_L_P_C_arm"]["seeds_LminusP"]["strata_signs"]},
       "interpretation":"linear lamp arrays vs point lights differ also in lamp number/footprint, not a geometry-only causal contrast",
       "original_lit_dark_joint_variance_on_randomized_field_scale":None},
      "Grognuz2026":{
       "source_type":"supplementary Table S3 derived block means (NOT original 425 per-pot observed seed rows)",
       "daytime_90min_pollen_access_only_after_night_bagging":True,
       "Trifolium":{"choice_blocks":t["n_blocks"],"dark_zero_blocks":len(t["dark_zero_blocks"]),
                    "positive":t["all_block_signs"]["lit_positive"],
                    "negative":t["all_block_signs"]["lit_negative"],
                    "mean_lit_minus_dark_seeds":t["mean_lit_minus_dark_seed_per_target"],
                    "excluding_dark_zero_posthoc":
                      t["dark_nonzero_only_posthoc"]},
       "Centaurea":{"choice_blocks":c["n_blocks"],"dark_zero_blocks":len(c["dark_zero_blocks"]),
                    "positive":c["all_block_signs"]["lit_positive"],
                    "negative":c["all_block_signs"]["lit_negative"],
                    "mean_lit_minus_dark_seeds":c["mean_lit_minus_dark_seed_per_target"],
                    "excluding_dark_zero_posthoc":
                      c["dark_nonzero_only_posthoc"]},
       "original_independent_lighting_sites_recovered":False},
      "mathematical_control_outcome_coupling":{
          "identity":"Cov(C, L-C) = Cov(C,L) - Var(C)",
          "under_independent_equal_null":"Cov(C,L-C) = -Var(C)",
          "strict_dark_zero_and_nonnegative_seed_count":"L-C >= 0 necessarily",
          "not_evidence_for_ALAN_rescue":True,
          "condition_on_dark_nonzero_is_post_treatment_observed_selection":True},
      "meta_inference":{
          "ALAN_overall_open_seed_lnRR":None,
          "ALAN_x_local_dark_night_dependence_V":None,
          "island_x_ALAN_reproduction_effect":None,
          "not_a_pooled_meta_analysis":True,
          "cross_study_signs_cannot_be_pooled_across_exposure_access_and_seed_denominators":True,
          "true_independent_lighting_sites_not_11_plus_28":True,
          "reasons":["different species, lighting regimes and original seed opportunity denominators",
                     "daytime-only after night bagging in Grognuz vs unbagged night exposure in Li",
                     "Li uses 3 true fields; Grognuz choice blocks share installations and pot reuse",
                     "selection on observed zero dark outcomes forces nonnegative comparisons",
                     "insufficient source covariance and same-population V for meta regression"]}}
    return summary

def main():
    liout=li.analyze(li.fetch())
    grout=grog.analyze(grog.load())
    reg=json.loads(SOURCE_REGISTER.read_text(encoding="utf-8"))
    out=result(liout,grout,reg)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_TWO_STUDY_FLOOR_GATE="+json.dumps(out,ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
