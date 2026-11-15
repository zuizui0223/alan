#!/usr/bin/env python3
"""Same-source daytime-only vs open reproduction: two distinct outcome stages.

Source is the historical alan night_exclusion/experiments.csv 42-row
curated literature table.  A matched study/site/year/plant/cultivar has
both a 'fruit set' and a 'seed set' endpoint, under the same access
contrast.  They are NOT two statistically independent effect studies.

Strict stage sensitivity omits:
  - fruit-set source values >1 (not a verified fruit/probability);
  - SmallFruit_2025 seed-quality rows, which do not match original fruit
    success counts in 9 of 20 group cells in the authentic Figshare;
  - explicit source unit-invalid rows.

Even for retained strata, seed-per-fruit can be conditional on forming
a fruit and the seed denominator differs from fruit probability.  Do not
make a direct causal net-seed-per-original-flower claim from the two ratios.
"""
from __future__ import annotations
import collections
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/"results/night_exclusion/experiments.csv"
OUTPUT=ROOT/"results/night_exclusion/stage_marginal_diel_pairs_20261115.json"
OUTPUT_CSV=ROOT/"results/night_exclusion/stage_marginal_strict_pairs_20261115.csv"
JOIN_KEYS=("plant","study_ID","locality","Site_period","year_start","plant_cultivar","additional_treatment")
TREAT=("fruit set","seed set")

def normalize(x:str)->str:
    return " ".join(str(x or "").strip().split())

def parse_rows(path=SRC):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def match_and_gate(rows):
    index=collections.defaultdict(dict)
    for row in rows:
        name=row["treatment_effectiveness_metric"]
        if name not in TREAT:continue
        key=tuple(normalize(row[c]) for c in JOIN_KEYS)
        if name in index[key]:
            raise AssertionError("Duplicate source outcome with same site-year-plant grouping")
        index[key][name]=row
    all_pairs=[]
    strict=[]
    rejected=[]
    for key,pair in sorted(index.items()):
        if set(pair)!=set(TREAT):continue
        fruit=pair["fruit set"];seed=pair["seed set"]
        vf=float(fruit["loss"]);vs=float(seed["loss"])
        fday=float(fruit["day_value"]);fopen=float(fruit["open_value"])
        sday=float(seed["day_value"]);sopen=float(seed["open_value"])
        if min(fday,fopen,sday,sopen)<0 or min(fopen,sopen)<=0:
            raise AssertionError("Original source measures cannot form these ratios")
        if not math.isclose(vf,1-fday/fopen,abs_tol=1e-9):
            raise AssertionError("Fruit source ratio drift")
        if not math.isclose(vs,1-sday/sopen,abs_tol=1e-9):
            raise AssertionError("Seed source ratio drift")
        p={"plant":key[0],"study_ID":key[1],"locality":key[2],
           "site_period":key[3],"year":key[4],"cultivar":key[5],
           "extra_treatment":key[6],
           "V_fruit_outcome":vf,"V_seed_outcome":vs,
           "V_seed_minus_V_fruit":vs-vf,
           "outcomes_nondisjoint_and_not_independent":True,
           "source_same_study_not_same_individual_plant_proven":True,
           "cluster_adjusted_SE":None}
        all_pairs.append(p)
        reason=None
        if key[1]=="SmallFruit_2025":
            reason="FRUIT_QUALITY_TABLE_ORIGINAL_9_OF_20_CELL_COUNT_MISMATCH"
        elif fday>1 or fopen>1:
            reason="FRUIT_SET_NOT_CERTIFIED_FRACTION_BECAUSE_ARM_VALUE_GT_ONE"
        elif "not used" in str(fruit["basis"]).lower() or "not used" in str(seed["basis"]).lower():
            reason="SOURCE_EXPLICITLY_INVALID_UNIT"
        if reason is None:strict.append(p)
        else:rejected.append({"plant":key[0],"study_ID":key[1],
                              "year":key[4],"reason":reason})
    if len(rows)==42 and (len(all_pairs),len(strict))!=(14,10):
        raise AssertionError("Review stage matching after source table version change")
    if len(rows)==42 and len(set(p["study_ID"] for p in strict))!=6:
        raise AssertionError("Unexpected strict source study count")
    return {
      "schema":"ALAN_B_DvN_stage_matched_outcome_gate_v1",
      "source":"results/night_exclusion/experiments.csv",
      "original_source_rows":len(rows),
      "all_stage_pairs":len(all_pairs),
      "strict_stage_pairs":len(strict),
      "strict_original_studies":len(set(p["study_ID"] for p in strict)),
      "rejected_pair_count":len(rejected),
      "rejected_reasons":rejected,
      "all_pair_seed_marginal_larger_count":sum(p["V_seed_minus_V_fruit"]>0 for p in all_pairs),
      "strict_seed_marginal_larger_count":sum(p["V_seed_minus_V_fruit"]>0 for p in strict),
      "all_pairs_not_independent_inference_units":True,
      "same_initial_flower_total_seed_fitness_identified":False,
      "selection_and_conditioning_on_formed_fruit_possible":True,
      "source_plant_random_effect_and_covariance_available":False,
      "source_species_guild_ALAN_impact_estimated":False,
      "geographic_island_ALAN_interaction_estimated":False,
      "global_meta_beta_V":None,
      "strict":strict,
      "all":all_pairs
    }

def main():
    out=match_and_gate(parse_rows())
    OUTPUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with OUTPUT_CSV.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(out["strict"][0]))
        writer.writeheader();writer.writerows(out["strict"])
    study_vals=collections.defaultdict(list)
    for p in out["strict"]:
        study_vals[p["study_ID"]].append(p["V_seed_minus_V_fruit"])
    study_means={s:sum(x)/len(x) for s,x in sorted(study_vals.items())}
    print(json.dumps({
      "original_rows":out["original_source_rows"],
      "paired_stage_comparisons":out["all_stage_pairs"],
      "source_quality_strict_pairs":out["strict_stage_pairs"],
      "source_quality_strict_studies":out["strict_original_studies"],
      "strict_seed_V_gt_fruit_V":out["strict_seed_marginal_larger_count"],
      "study_level_mean_stage_differences":study_means,
      "global_meta_beta_V":None
    },ensure_ascii=False))

if __name__=="__main__":
    main()
