#!/usr/bin/env python3
"""Fail closed if original EERES raw XLSX gives changed 2x2 ALAN×flow evidence."""
import argparse
import json
import math
from pathlib import Path

def verify(frozen,original):
    root=original.get("eeres",{})
    x=root.get("flow_moderator_exploratory",{})
    if original.get("status")!="SOURCE_AUTHENTICATED_EXPLORATORY_TWO_PROGRAMME_COMPARISON_NOT_ACCEPTED_META":
        raise ValueError("STOP: no genuine cross-programme source replay")
    if root.get("source_sha256")!=frozen["source_original_workbook_sha256"]:
        raise ValueError("STOP: changed original author workbook SHA")
    if root.get("physical_header_n")!=6 or x.get("n_true_independent_header_tanks")!=6:
        raise ValueError("STOP: missing six physical independent blocks")
    checks={
        "normal_flow_lit_vs_dark_RRR":root["four_order_composition_RRR"],
        "reduced_flow_lit_vs_dark_RRR":x["reduced_flow_light_composition_RRR"],
        "reduced_flow_lnRRR":x["reduced_flow_logRRR"],
        "reduced_flow_lnRRR_se":x["reduced_flow_delta_SE"],
        "source_estimate_moderation_lnR_of_RRRs":x["within_same_six_HT_moderation_logRRR"],
        "source_estimate_moderation_R_of_RRRs":x["moderation_ratio_of_RRRs"],
        "source_estimate_moderation_SE_full_shared_header":x["joint_shared_HT_delta_SE"],
    }
    for name,val in checks.items():
        if not math.isclose(float(val),float(frozen[name]),rel_tol=1e-11,abs_tol=1e-12):
            raise ValueError(f"STOP: original source mismatch: {name}")
    for f,orig in (
        ("reduced_flow_descriptive_t5_95",x["reduced_flow_descriptive_t5_95"]),
        ("source_estimate_moderation_descriptive_t5_95",x["descriptive_t5_95_moderator"])
    ):
        if len(orig)!=2 or any(not math.isclose(float(a),float(b),rel_tol=1e-11,abs_tol=1e-12)
                               for a,b in zip(frozen[f],orig)):
            raise ValueError(f"STOP: source interval mismatch: {f}")
    if original.get("meta_k_added")!=0 or not x.get("not_a_water_wet_vs_dry_gradient_or_plant_fitness"):
        raise ValueError("STOP: new meta k or unsupported plant mechanism")
    if not math.isclose(float(frozen["source_estimate_moderation_R_of_RRRs"]),
                        float(frozen["reduced_flow_lit_vs_dark_RRR"])/
                        float(frozen["normal_flow_lit_vs_dark_RRR"]),rel_tol=1e-12):
        raise ValueError("STOP: factorial moderator RR inconsistent")
    return {"status":"ORIGINAL_SOURCE_REPLAY_MATCHES_FROZEN_FLOW_MODERATOR",
            "n_physical_HT":6,"light_x_flow_interaction_RRR":frozen["source_estimate_moderation_R_of_RRRs"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--frozen",type=Path,required=True)
    ap.add_argument("--source-replay",type=Path,required=True)
    args=ap.parse_args()
    a=json.loads(args.frozen.read_text(encoding="utf-8"))
    b=json.loads(args.source_replay.read_text(encoding="utf-8"))
    print(json.dumps(verify(a,b),ensure_ascii=False))
if __name__=="__main__":
    main()
