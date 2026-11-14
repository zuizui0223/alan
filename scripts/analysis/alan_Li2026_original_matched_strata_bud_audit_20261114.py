#!/usr/bin/env python3
"""Li (2026) authentic-source matched stratum and starting-bud audit.

Disentangle original three-site spatial replication, species and deployment-round
opportunities in Li et al. DOI 10.64898/2026.07.10.737511. This cannot
identify pure lamp GEOMETRY because the author array has 26 bulbs, point has 4.
No plant or site raw rows are exported; source Github blob SHA must match.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from scripts.analysis.alan_Li2026_real_geometry_seed_hayfield_source_20261114 import (
    source_bytes, parse, summary, contrast, require, GIT_BLOB_SHA1
)

ARMS=("L","P","C")
PAIR=(("L","C"),("L","P"),("P","C"))
SPP=("ML","SN","ZC")

def valid(rows):
    return [r for r in rows if r["pollinated"]!="NA"]

def eligible_group(rows, arm):
    v=valid([r for r in rows if r["Treatment"]==arm])
    if not v:
        return None
    require(all(r["Starting"] not in ("NA","") for r in v),
            "Baseline starting buds unavailable")
    require(all(int(r["Starting"])>=0 for r in v),
            "Starting bud counts changed")
    n=len(v)
    s=sum(int(r["Starting"]) for r in v)
    return {"n":n,"seed_positive":sum(int(r["pollinated"]) for r in v),
        "seed_presence_fraction":sum(int(r["pollinated"]) for r in v)/n,
        "sum_starting_buds":s, "mean_starting_buds_per_plant":s/n,
        "seeds_per_starting_bud":sum(int(r["Seeds"]) for r in v)/s if s>0 else None,
        "fruits_per_starting_bud":sum(int(r["Fruit"]) for r in v)/s if s>0 else None,
        "sum_seeds":sum(int(r["Seeds"]) for r in v),
        "sum_fruits":sum(int(r["Fruit"]) for r in v),
    }

def pair_result(rows, treatment, control):
    a=eligible_group(rows,treatment); b=eligible_group(rows,control)
    if a is None or b is None:
        return {"comparable":False,
                "missing_arms":[x for x,y in ((treatment,a),(control,b)) if y is None]}
    return {"comparable":True,"treatment":treatment,"control":control,
        "seed_presence_difference":a["seed_presence_fraction"]-b["seed_presence_fraction"],
        "starting_bud_mean_difference":a["mean_starting_buds_per_plant"]-b["mean_starting_buds_per_plant"],
        "seeds_per_starting_bud_difference":
           a["seeds_per_starting_bud"]-b["seeds_per_starting_bud"] if a["seeds_per_starting_bud"] is not None and b["seeds_per_starting_bud"] is not None else None,
        "fruits_per_starting_bud_difference":
           a["fruits_per_starting_bud"]-b["fruits_per_starting_bud"] if a["fruits_per_starting_bud"] is not None and b["fruits_per_starting_bud"] is not None else None,
        "arm_summaries":{treatment:a,control:b}
    }

def strata_table(rows, groupfields, treatment, control):
    keys=sorted({tuple(r[g] for g in groupfields) for r in rows})
    items=[]
    for key in keys:
        sub=[r for r in rows if tuple(r[g] for g in groupfields)==key]
        item=pair_result(sub,treatment,control)
        item["stratum"]=dict(zip(groupfields,key))
        items.append(item)
    supported=[x for x in items if x["comparable"]]
    d=[x["seed_presence_difference"] for x in supported]
    r=[x["seeds_per_starting_bud_difference"] for x in supported
       if x["seeds_per_starting_bud_difference"] is not None]
    return {"group_fields":list(groupfields),"pair":treatment+"-"+control,
        "strata":items, "total_observed_strata":len(items),
        "n_comparable_strata":len(supported),"n_unmatched":len(items)-len(supported),
        "n_positive":sum(x>0 for x in d),
        "n_negative":sum(x<0 for x in d),
        "n_equal":sum(x==0 for x in d),
        "equal_supported_stratum_mean_risk_difference":sum(d)/len(d) if d else None,
        "equal_supported_stratum_mean_seeds_per_starting_bud_difference":
           sum(r)/len(r) if r else None,
        "not_independent_experimental_replicates":True}

def leave_species_out(rows, treat, control):
    out=[]
    for removed in SPP:
        remain=[r for r in rows if r["species"]!=removed]
        out.append({"omitted_species":removed,
                    "pooled_plant_difference":pair_result(remain,treat,control)})
    return out

def by_round(rows,treat,control):
    # Round 1 is only a within-time source cohort; ML was deployed in all
    # three rounds, SN/ZC in the first only, so round differences are
    # confounded with species composition.
    return [{"round":q,"same_round_original_plant_level":pair_result(
        [r for r in rows if int(r["round"])==q],treat,control)}
        for q in (1,2,3)]

def source_coverage(rows):
    by=Counter(r["species"] for r in rows)
    # Author preprint METHODS states potential deployments ML 54, SN 18, ZC 18.
    # Published CSV rows not identical to potential deployment slots.
    intended={"ML":54,"SN":18,"ZC":18}
    return {"author_methods_potential_plant_deployment_slots":
        intended,"total_potential_slots_described":sum(intended.values()),
        "source_data_rows_by_species":dict(by),
        "total_original_csv_rows":len(rows),
        "difference_potential_minus_archived_source_rows":
            {s:intended[s]-by.get(s,0) for s in intended},
        "do_not_infer_unrecorded_rows_as_failed_pollination":True,
        "no_source_verified_individual_missingness_reason":True}

def missingness_tipping_bound(rows, treatment, control, max_missing=40):
    """Extreme unfavorable hypothetical allocation of unobserved plant slots.

    For each missing M, allocate x failures to treatment and M-x successes
    to control. This ignores species, site, round and actual eligibility;
    it is only a mathematical sensitivity, never source-attributed outcomes.
    """
    a=eligible_group(rows,treatment)
    b=eligible_group(rows,control)
    require(a is not None and b is not None,"Both source groups required")
    a_yes=a["seed_positive"];a_n=a["n"]
    b_yes=b["seed_positive"];b_n=b["n"]
    def find(m):
        options=[]
        for x in range(m+1):
            y=m-x
            delta=a_yes/(a_n+x)-(b_yes+y)/(b_n+y)
            options.append((delta,x,y))
        diff,x,y=min(options,key=lambda q:q[0])
        return {"hypothetical_total_missing":m,
                "minimum_possible_unadjusted_risk_difference":diff,
                "unobserved_treatment_failures":x,
                "unobserved_control_successes":y}
    original_delta=a_yes/a_n-b_yes/b_n
    bound11=find(11)
    reversal=next((find(m) for m in range(max_missing+1)
                   if find(m)["minimum_possible_unadjusted_risk_difference"]<=0),None)
    return {"pair":treatment+"-"+control,
            "original_raw_difference":original_delta,
            "hypothetical_11_missing_worst_case":bound11,
            "first_possible_sign_reversal_if_adversarial_missing":reversal,
            "not_actual_assignment_or_missingness_model":True,
            "treatment_specific_missing_counts_unknown":True,
            "valid_only_as_unstratified_algebra_not_per_field_causal_identification":True}

def analyze(rows, raw):
    results={"schema":"ALAN_Li2026_79_original_plants_matched_species_field_round_baseline_bias_v1",
      "source":{"original_author_csv_sha256":hashlib.sha256(raw).hexdigest(),
                "original_author_git_blob_sha1":GIT_BLOB_SHA1,
                "study_DOI":"10.64898/2026.07.10.737511",
                "original_independent_field_blocks":3,
                "experimental_treatment_bulbs_by_author_method":{"L":26,"P":4,"C":0},
                "L_and_P_photometric_power_and_exposure_geometry_not_held_fixed":True,
                "source_preprint_not_peer_reviewed":True,
                "published_methods_moth_visitors_not_directly_observed":True},
      "source_data_vs_methods":source_coverage(rows),
      "starting_bud_baselines":{
        "whole_source_by_arm":{arm:eligible_group(rows,arm) for arm in ARMS},
        "original_author_methods_report_no_significant_baseline_bud_treatment_difference":True,
        "do_not_infer_starting_bud_difference_statistically_significant":True},
      "raw_pooled_three_arm_pairs":{a+"-"+b:pair_result(rows,a,b) for a,b in PAIR},
      "adversarial_source_data_vs_90_potential_slots_bound":{
        a+"-"+b:missingness_tipping_bound(rows,a,b) for a,b in PAIR},
      "matched_original_species_x_field":{
        a+"-"+b:strata_table(rows,("species","replicate"),a,b) for a,b in PAIR},
      "matched_original_species_x_field_x_round":{
        a+"-"+b:strata_table(rows,("species","replicate","round"),a,b) for a,b in PAIR},
      "source_same_round_responses":{
        a+"-"+b:by_round(rows,a,b) for a,b in PAIR},
      "whole_source_leave_one_species_out":{
        a+"-"+b:leave_species_out(rows,a,b) for a,b in PAIR},
      "limits":[
       "The 79 available plant rows are not 79 independent lighting-treatment installations",
       "Original source paper describes 90 potential deployments, but raw source has 79 rows: the 11-row difference cannot be turned into inferred zero outcomes",
       "Array has 26 bulbs, point light has 4, control zero, exposure intensity/footprint and geometry confounded",
       "Original species and round strongly confounded: two species only sampled round 1",
       "Starting buds are a pre-exposure opportunity proxy, NOT a measured local flower opening or pollen deposition",
       "Ratio of total seeds to initial buds is a sensitivity measure, NOT per-ovule true fertilization success",
       "Any within-field/species strata are repeated descriptions of only three independent fields",
       "Author original GLM did not detect array-vs-point pairwise difference: do not infer a significant geometry effect from a subgroup",
       "No nocturnal moth visitation, effective stigma pollen transfer, light-corridor movement, or island comparison measured"
      ]
    }
    # Source sign/invariance gates, not desired positive biological claims.
    require(source_coverage(rows)["total_original_csv_rows"]==79,"Original CSV size changed")
    baseline=results["starting_bud_baselines"]["whole_source_by_arm"]
    require(baseline["L"]["n"]==27 and baseline["C"]["n"]==26 and baseline["P"]["n"]==25,"Source arm counts changed")
    raw_lc=results["raw_pooled_three_arm_pairs"]["L-C"]
    require(abs(raw_lc["seed_presence_difference"]-(20/27-8/26))<1e-12,"Original source outcome changed")
    strata=results["matched_original_species_x_field"]["L-P"]
    require(strata["n_comparable_strata"]==9 and
            strata["n_positive"]==4 and strata["n_negative"]==1 and
            strata["n_equal"]==4,"Original L-v-P within-species field discordance changed")
    return results

def main(out):
    out=Path(out);require(not out.exists(),"Refuse overwrite")
    raw=source_bytes(); rows=parse(raw)
    r=analyze(rows,raw)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_Li2026_MATCHED_STRATA_ORIGINAL",json.dumps(r,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
    main(p.parse_args().out)
