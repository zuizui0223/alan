#!/usr/bin/env python3
"""Species-level MOIST robustness on the immutable 3,901-species alan/main source.

This is *selection/trait-coverage diagnostic* of a previously selected cohort.
It is NOT an all-species continuous ALAN-response refit and does not identify a causal
wetland, pollinator, predation or flowering pathway.
"""
import csv
import gzip
import hashlib
import json
import math
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "results/all_species_forest/decline_traits"
LOCKED_MAIN_SHA = "52ff636fa8ffc56e08fce5f7a47795af90099da9"
SEED = 20261117

def read_source():
    path = SOURCES / "moisture_evidence.csv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        lines = list(csv.DictReader(handle))
    assert len(lines) == 3901, f"expected 3,901 species, observed {len(lines)}"
    assert len({r["sp"] for r in lines}) == 3901, "duplicate original species key"
    assert sum(r["candidate"] == "True" for r in lines) == 283
    assert all(r["candidate"] in ("True", "False") for r in lines)
    assert all(set(r) >= {"sp", "species", "candidate", "indicator",
                          "site_region", "island_band", "site_band"} for r in lines)
    return lines, hashlib.sha256(path.read_bytes()).hexdigest()

def get_strata(data):
    groups = defaultdict(lambda: {"case": [], "ref": [], "case_all": 0, "ref_all": 0})
    for r in data:
        groupkey = (r["site_region"], r["island_band"], r["site_band"])
        if any(v in ("", "nan") for v in groupkey):
            continue
        g = groups[groupkey]
        label = r["candidate"] == "True"
        g["case_all" if label else "ref_all"] += 1
        if r["indicator"] not in ("", "nan", "NaN"):
            val = float(r["indicator"])
            if not math.isfinite(val):
                raise AssertionError("nonfinite EIV trait")
            g["case" if label else "ref"].append(val)
    return {key: g for key, g in groups.items() if g["case"] and g["ref"]}, groups

def difference(groups):
    n, signed = 0, 0.
    for g in groups.values():
        a,b = g["case"], g["ref"]
        n += len(a)
        signed += sum(a)-len(a)*sum(b)/len(b)
    return (signed/n if n else None), n

def category(val):
    return "low_le7" if val <= 7 else ("mid_7to9" if val <= 9 else "verywet_gt9")

def shape(groups):
    total = sum(len(g["case"]) for g in groups.values())
    observed = {cat: 0. for cat in ["low_le7","mid_7to9","verywet_gt9"]}
    expected = dict.fromkeys(observed, 0.)
    mh_num = mh_den = 0.
    for g in groups.values():
        a,b = g["case"],g["ref"]
        na,nb=len(a),len(b)
        for v in a: observed[category(v)] += 1.
        for cat in expected:
            expected[cat] += na * sum(category(v)==cat for v in b)/nb
        ca = sum(v>7 for v in a)
        cr = sum(v>7 for v in b)
        n = na + nb
        mh_num += (ca * (nb-cr)) / n
        mh_den += ((na-ca) * cr) / n
    assert total > 0 and mh_num > 0 and mh_den > 0
    return {"matched_case_count": total,
            "bins": {cat:{"candidate_n":observed[cat],
                          "matched_ref_expected_n":expected[cat],
                          "candidate_percent":100*observed[cat]/total,
                          "matched_reference_percent":100*expected[cat]/total,
                          "percentage_point_difference":100*(observed[cat]-expected[cat])/total}
                     for cat in observed},
            "wet_gt7_conditional_MH_odds_ratio":mh_num/mh_den}

def missingness(data, matched, full):
    total_n=len(data)
    case=[r for r in data if r["candidate"]=="True"]
    ref=[r for r in data if r["candidate"]=="False"]
    known=lambda s:sum(r["indicator"] not in ("","nan","NaN") for r in s)
    ret={"all_species":total_n,"candidate_count":len(case),"reference_count":len(ref),
         "candidate_indicator_known":known(case),"reference_indicator_known":known(ref),
         "candidate_coverage":known(case)/len(case),"reference_coverage":known(ref)/len(ref),
         "matched_strata_count":len(matched)}
    # Count EIV-unknown species in the *same* strata admitted by the existing
    # positive complete-case analysis. Missing labels remain missing, not zero.
    ret.update({
        "matched_strata_case_all":sum(g["case_all"] for g in matched.values()),
        "matched_strata_ref_all":sum(g["ref_all"] for g in matched.values()),
        "matched_strata_case_known":sum(len(g["case"]) for g in matched.values()),
        "matched_strata_ref_known":sum(len(g["ref"]) for g in matched.values())
    })
    ret["case_coverage_in_matched_strata"] = (
        ret["matched_strata_case_known"]/ret["matched_strata_case_all"])
    ret["reference_coverage_in_matched_strata"] = (
        ret["matched_strata_ref_known"]/ret["matched_strata_ref_all"])
    return ret

def geography(groups):
    by_region=defaultdict(dict)
    for key,g in groups.items():by_region[key[0]][key]=g
    point,n=difference(groups)
    stats=[]
    for region,subset in by_region.items():
        d,k=difference(subset)
        stats.append({"region":region,"candidates":k,"difference":d})
    stats.sort(key=lambda r:(-r["candidates"],r["region"]))
    assert sum(r["candidates"] for r in stats)==n
    leave_one=[]
    for x in stats:
        d,k=difference({key:g for key,g in groups.items() if key[0]!=x["region"]})
        leave_one.append({"omitted":x["region"],"remaining_cases":k,"difference":d})
    omit={stats[0]["region"], stats[1]["region"]}
    dbl,k=difference({key:g for key,g in groups.items() if key[0] not in omit})
    region_equal=sum(r["difference"] for r in stats)/len(stats)
    # A regional sign test is DESCRIPTIVE because the 12 sectors are not independent
    # experimental interventions and cases were preselected from these time series.
    positives=sum(r["difference"]>0 for r in stats)
    # P(Binomial(12,.5) is at least as extreme as positives); descriptive only.
    pos_tail=sum(math.comb(len(stats),i) for i in range(positives,len(stats)+1))
    neg_tail=sum(math.comb(len(stats),i) for i in range(0,len(stats)-positives+1))
    p_sign=min(1.,2*min(pos_tail,neg_tail)/2**len(stats))
    rng=random.Random(SEED)
    # Resample 12 entire geographic origin groups. This retains within-origin
    # rows and is not a correction for adaptive candidate selection.
    region_stats={r["region"]:(r["difference"],r["candidates"]) for r in stats}
    rkeys=list(region_stats)
    draws=[]
    for _ in range(3999):
        sample=[region_stats[rng.choice(rkeys)] for _ in rkeys]
        draws.append(sum(diff*n for diff,n in sample)/sum(n for diff,n in sample))
    draws.sort()
    return {"n_regions":len(stats),"positive":positives,"negative":len(stats)-positives,
            "region_equal_mean_difference":region_equal,
            "region_sign_test_descriptive_two_sided_p":p_sign,
            "largest_two_region_ids":[r["region"] for r in stats[:2]],
            "largest_two_case_count":sum(r["candidates"] for r in stats[:2]),
            "excluding_largest_two_difference":dbl,"excluding_largest_two_cases":k,
            "single_region_leave_out_min":min(r["difference"] for r in leave_one),
            "single_region_leave_out_max":max(r["difference"] for r in leave_one),
            "12_region_resampling_descriptive_interval":[draws[99],draws[3899]],
            "by_region":stats,"leave_one":leave_one}

def run():
    src,sha=read_source()
    matched,full=get_strata(src)
    d,n=difference(matched)
    assert len(matched)==49, f"matched strata changed: {len(matched)}"
    assert n==158 and abs(d-0.5922104455502863)<1e-9, (n,d)
    return {
        "schema":"alan_wet_moisture_species_support_v1",
        "main_source_commit":LOCKED_MAIN_SHA,
        "source_file":"results/all_species_forest/decline_traits/moisture_evidence.csv.gz",
        "source_gzip_sha256":sha,
        "interpretation":"previously selected candidates; descriptive within-stratum species trait audit, NOT full 3901 continuous light response",
        "frozen_primary_moisture_difference":d,
        "coverage":missingness(src,matched,full),
        "moisture_bins_and_odds":shape(matched),
        "geographic_sensitivity":geography(matched),
        "not_estimated":["continuous_3901_species_ALAN_slope","wetness_x_ALAN_causal_effect",
                         "independent_reproduction_meta_effect"]
    }

if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2,sort_keys=True))
