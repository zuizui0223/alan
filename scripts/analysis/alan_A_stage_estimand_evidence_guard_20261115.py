#!/usr/bin/env python3
"""Validate source-defined ALAN reproduction stage identities; NEVER pool scales."""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
MAP=ROOT/"config/alan_A_stage_estimand_evidence_map_20261115_v2.json"
REG=ROOT/"config/alan_A_reproductive_screen_20261115.json"
def validate(stage_map,source_reg):
    papers=stage_map["papers"]
    reg={r["id"]:r for r in source_reg["entries"]}
    assert len(papers)==stage_map["paper_count_expected"]==9
    assert len({p["id"] for p in papers})==9
    assert {p["id"] for p in papers}=={
        r["id"] for r in source_reg["entries"] if
        r["a_primary_stage"]!="exclude_A_no_reproduction"
        and r["light_design"].startswith("manipulated")
    }, "Do not silently add observational/phenology or duplicated papers"
    metrics=Counter()
    for p in papers:
        assert p["doi"]==reg[p["id"]]["doi"]
        assert p["outcomes"]
        for x in p["outcomes"]:
            assert x["scope"] and x["contrast"] and x["scale"] and x["stage"]
            assert x["ready"] is False, "No cluster-harmonized joint variance registered"
            assert x["ci95"] is None or len(x["ci95"])==2
            if x["scale"] in ("direction_only","not_reported"):
                assert x["value"] is None
            if x["scale"] in ("published_RR","descriptive_RR",
                             "derived_ratio_from_published_stage_means",
                             "original_aggregate_RR","RR_approx","OR_conditional"):
                assert x["value"] is not None and x["value"]>0
            if x["scale"]=="author_model_beta_standardized":
                assert x.get("ci85") is not None and x.get("ci95") is None
            metrics[x["stage"]]+=1
    assert len(stage_map["additional_grey_leads"])==1
    gray=stage_map["additional_grey_leads"][0]
    assert gray["status"]=="conference_abstract_only_not_meta_ready"
    assert gray["exact_effect"] is None and gray["independent_lamp_units"] is None
    # Thesis/preprint/abstract and multiple endpoint rows must not increase k.
    assert not any(x["ready"] for p in papers for x in p["outcomes"])
    assert next(p for p in papers if p["id"]=="lu2026jpe")["outcomes"][1]["value"] is None
    assert next(p for p in papers if p["id"]=="macgregor2019")["outcomes"][0]["value"] is None
    return {"direct_manipulated_papers":len(papers),
            "source_effect_rows_not_independent":sum(len(p["outcomes"]) for p in papers),
            "different_stages":len(metrics),
            "confirmed_ALAN_x_independent_dark_V_meta_effects":0,
            "additional_unverified_conference_leads":len(stage_map["additional_grey_leads"]),
            "scale_categories":dict(sorted(Counter(x["scale"] for p in papers for x in p["outcomes"]).items())),
            "stage_categories":dict(sorted(metrics.items())),
            "is_meta_analysis":False}
def main():
    a=json.loads(MAP.read_text(encoding="utf-8"))
    b=json.loads(REG.read_text(encoding="utf-8"))
    print(json.dumps(validate(a,b),ensure_ascii=False,indent=2))
if __name__=="__main__":main()
