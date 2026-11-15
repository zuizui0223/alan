#!/usr/bin/env python3
"""Original French diel fruit-set arm contrasts with conservative small-n guards.

Source: Blareau & Requier 2025, DOI 10.1111/jen.70019; Figshare
DOI 10.6084/m9.figshare.29930195.v1.

Each of the original 80 fruit-set CSV rows is an INFLORESCENCE with successful
and unsuccessful fruits. Inflorescences are NOT independent plants. Plant IDs
are absent, so any bootstrap here is descriptive/anti-conservative about true
between-plant uncertainty. Main B result is a point estimate per taxon-year;
only 2024 arms >=3 inflorescences have a labelled exploratory interval.
"""
from __future__ import annotations
import collections
import importlib.util
import json
import pathlib
import random

BASE = pathlib.Path(__file__).resolve().parent
S = importlib.util.spec_from_file_location("alan_crop_original",BASE/"alan_P1_original_smallfruit_diel_dependence_20261115.py")
ORIG = importlib.util.module_from_spec(S)
S.loader.exec_module(ORIG)
SEED = 20261115
BOOT = 20000
MAP = {"BC":"Ribes nigrum","RB":"Rubus idaeus","RD":"Ribes rubrum"}
ARMS = ("closed","day_only","night_only","day_night_open")


def grouped_inflorescences(rows: list[dict]) -> dict[tuple[str,str,str],list[tuple[int,int]]]:
    out=collections.defaultdict(list)
    for row in rows:
        key=(str(row["Plant_sp"]).strip(),str(row["Year"]).strip(),ORIG.normalize_treatment(row["Treatment"]))
        if key[0] not in MAP or key[2] not in ARMS:
            raise ValueError("Unexpected source plant or treatment: "+repr(key))
        success=ORIG.int_count(row["Fruit_set_success"])
        failed=ORIG.int_count(row["Fruit_set_failure"])
        if success+failed <1:
            raise ValueError("Empty original inflorescence")
        out[key].append((success,failed))
    return dict(out)


def aggregate_probability(inflorescences: list[tuple[int,int]]) -> float | None:
    successes=sum(x[0] for x in inflorescences)
    total=sum(sum(x) for x in inflorescences)
    return None if not total else successes/total


def dependence(day: list[tuple[int,int]], all_access: list[tuple[int,int]]) -> float | None:
    pd=aggregate_probability(day);po=aggregate_probability(all_access)
    return None if po is None or po==0 or pd is None else 1-pd/po


def pct(sorted_samples: list[float], p: float) -> float:
    k=(len(sorted_samples)-1)*p
    lo=int(k);hi=min(lo+1,len(sorted_samples)-1)
    return sorted_samples[lo]*(hi-k)+sorted_samples[hi]*(k-lo)


def diagnostics(day: list[tuple[int,int]],all_access: list[tuple[int,int]],seed: int, boot_n: int=BOOT) -> dict:
    point=dependence(day,all_access)
    if point is None:
        return {"V":None,"CI_descriptive":None,"reason":"OPEN_RATE_ZERO_OR_NO_VALID_PAIR"}
    n_day,n_open=len(day),len(all_access)
    out={"V":point,"n_day_inflorescences":n_day,"n_open_inflorescences":n_open,
         "bootstrap_unit":"recorded_inflorescence_not_independent_plant",
         "ci_status":"EXPLORATORY_ONLY_PLANT_IDS_MISSING",
         "CI_descriptive":None,"bootstrap_draws_valid":0,"bootstrap_draws_invalid":0,
         "leave_one_inflorescence_min":None,"leave_one_inflorescence_max":None,
         "leave_one_signs":None,
         "true_plant_cluster_CI":None}
    leave=[]
    if n_day>1:
        leave += [dependence(day[:i]+day[i+1:],all_access) for i in range(n_day)]
    if n_open>1:
        leave += [dependence(day,all_access[:i]+all_access[i+1:]) for i in range(n_open)]
    leave=[x for x in leave if x is not None]
    if leave:
        out["leave_one_inflorescence_min"]=min(leave)
        out["leave_one_inflorescence_max"]=max(leave)
        out["leave_one_signs"]={"positive":sum(x>0 for x in leave),
                                "negative":sum(x<0 for x in leave),
                                "zero":sum(x==0 for x in leave)}
    if n_day<3 or n_open<3:
        out["ci_status"]="SUPPRESSED_LT3_INFLORESCENCES_PER_ARM"
        return out
    rng=random.Random(seed)
    arr=[]
    for _ in range(boot_n):
        sample_day=[day[rng.randrange(n_day)] for _ in range(n_day)]
        sample_open=[all_access[rng.randrange(n_open)] for _ in range(n_open)]
        val=dependence(sample_day,sample_open)
        if val is None:
            out["bootstrap_draws_invalid"]+=1
        else:
            arr.append(val)
    out["bootstrap_draws_valid"]=len(arr)
    if len(arr)>=.9*boot_n:
        arr.sort()
        out["CI_descriptive"]=[pct(arr,.025),pct(arr,.975)]
        out["bootstrap_fraction_V_gt0"]=sum(x>0 for x in arr)/len(arr)
        out["bootstrap_median"]=pct(arr,.5)
        out["ci_status"]="NONPARAMETRIC_INFLORESCENCE_RESAMPLING_NOT_PLANT_LEVEL_CI"
    else:
        out["ci_status"]="SUPPRESSED_TOO_MANY_BOOTSTRAP_ZERO_OPEN_RATE"
    return out


def analyze(rows: list[dict], boot_n: int=BOOT) -> dict:
    groups=grouped_inflorescences(rows)
    combos=sorted(set((g[0],g[1]) for g in groups))
    if len(rows)==80 and len(combos)!=5:
        raise AssertionError("Expected 5 original species-year groups")
    all_out=[]
    for i,(taxon,year) in enumerate(combos):
        arms={arm:groups.get((taxon,year,arm),[]) for arm in ARMS}
        probs={arm:aggregate_probability(records) for arm,records in arms.items()}
        counts={arm:{"inflorescences":len(records),
                     "success":sum(x[0] for x in records),
                     "failure":sum(x[1] for x in records)}
                for arm,records in arms.items()}
        d=diagnostics(arms["day_only"],arms["day_night_open"],SEED+100*i,boot_n=boot_n)
        p_open=probs["day_night_open"]
        num=probs["closed"];night=probs["night_only"]
        V_night_only_relative_open=None if not p_open else night/p_open
        C_closed_relative_open=None if not p_open else num/p_open
        expected_independent=None
        if num is not None and num<1 and night is not None and probs["day_only"] is not None:
            expected_independent=1-(1-probs["day_only"])*(1-night)/(1-num)
        nonadditivity=None if expected_independent is None or p_open is None else p_open-expected_independent
        all_out.append({
          "taxon_code":taxon,"taxon":MAP[taxon],"year":int(year),
          "arm_success_fail_inflorescence_counts":counts,
          "arm_fruit_success_probability":probs,
          "night_access_marginal_V":d,
          "night_only_over_open":V_night_only_relative_open,
          "closed_over_open":C_closed_relative_open,
          "nocturnal_vs_day_night_difference":None if p_open is None or night is None else night-p_open,
          "independent_channel_open_probability_hypothetical":expected_independent,
          "open_minus_independent_hypothetical":nonadditivity,
          "caveat":"Independence baseline is illustrative; differences are NOT proof of pollinator facilitation, complementarity, compensatory response, or causal nonadditivity."
        })
    return {
      "schema":"ALAN_P1_B1_2025_French_smallfruit_sensitivity_v1",
      "source":"DOI 10.1111/jen.70019; original hash-pinned Figshare DOI 10.6084/m9.figshare.29930195.v1",
      "source_rows":len(rows),"primary_analysis":"descriptive marginal nocturnal access dependence in fruit formation",
      "sample_unit":"inflorescence, repeated within plants and sites; plants IDs absent",
      "site_year_structure":"2020 and 2024 different sites; red and black currants have different-year observations from related propagation material; years are not independent randomized ALAN experiments.",
      "no_ALAN_manipulated":True,"moth_pollinator_identification":False,
      "plant_random_effect_fit":False,"ALAN_meta_beta_V":None,
      "uncertainty":"Leave-one-inflorescence sensitivity all five; exploratory within-arm inflorescence bootstrap only for groups with >=3 inflorescences each. No valid plant-level 95% CI.",
      "seed_quality":"Separate original 223 formed-fruit rows cannot be joined 1:1 to all 258 fruit success counts in nine taxon-year-treatment cells. Do not compute total seed production per initial flower.",
      "boot_reps":boot_n,"seed":SEED,
      "results":all_out
    }


def main():
    orig=ORIG.read_original("fruit_set")
    out=analyze(orig)
    dest=pathlib.Path("results/alan_P1_B1_smallfruit_inflorescence_sensitivity_20261115.json")
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"source_rows":out["source_rows"],"boot_reps":out["boot_reps"],
        "study_year_results":[{"taxon":r["taxon"],"year":r["year"],
        "V":r["night_access_marginal_V"]["V"],
        "LOIO":[r["night_access_marginal_V"]["leave_one_inflorescence_min"],
                r["night_access_marginal_V"]["leave_one_inflorescence_max"]],
        "bootstrap_interval":r["night_access_marginal_V"]["CI_descriptive"],
        "ci_status":r["night_access_marginal_V"]["ci_status"],
        "open_minus_null":r["open_minus_independent_hypothetical"]}
            for r in out["results"]]},ensure_ascii=False))


if __name__=="__main__":
    main()
