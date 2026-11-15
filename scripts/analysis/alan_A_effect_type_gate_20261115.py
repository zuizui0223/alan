#!/usr/bin/env python3
"""ALAN A: explicitly incompatible source-reported reproductive effect measures.

Use source-published aggregate numbers only. This is NOT an inferential meta-fit.
Never mislabel a logistic odds ratio or conditional reproduction index as
a per-initial-flower log response ratio.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

def log_response_ratio(x_lit: float, x_dark: float) -> float:
    if not (math.isfinite(x_lit) and math.isfinite(x_dark) and x_lit > 0 and x_dark > 0):
        raise ValueError("Positive finite treatment means are required")
    return math.log(x_lit / x_dark)

def publication_compatible_examples() -> dict:
    li_p_lit, li_p_dark = 20/27, 8/26
    li_logrr = log_response_ratio(li_p_lit, li_p_dark)
    giavi_logit, giavi_se = -1.37, 0.726
    lu_seed_factor = 1 - 0.56
    out = {
        "schema": "ALAN_A_published_stage_compatibility_v1",
        "date": "2026-11-15",
        "not_meta_analysis": True,
        "confirmed_A_beta_V": None,
        "studies": [
            {
                "id": "li2026",
                "doi": "10.64898/2026.07.10.737511",
                "species": "3 species pooled over 3 fields",
                "outcome": "at least one seed per observed plant",
                "exposure": "linear light array versus no light",
                "n_success_lit": 20, "n_observed_lit": 27,
                "n_success_dark": 8, "n_observed_dark": 26,
                "light_dark_ratio": li_p_lit/li_p_dark,
                "ln_RR": li_logrr,
                "ln_RR_SE_independent_clustering": None,
                "primary_initial_flower_endpoint": False,
                "limitations": "Not initial flower fruit set; 3 shared experimental fields and 3 pooled species. No fabricated independent plant-level SE."
            },
            {
                "id": "lu2026jpe_pedicularis",
                "doi": "10.1093/jpe/rtag160",
                "species": "Pedicularis kansuensis",
                "outcome": "author-reported seed formation reduction, exact denominator and SE not yet recovered",
                "exposure": "ALAN versus natural dark",
                "light_dark_factor_from_author_percent": lu_seed_factor,
                "ln_ratio_of_author_percent": math.log(lu_seed_factor),
                "ln_RR_SE": None,
                "primary_initial_flower_endpoint": None,
                "limitations": "The published 56% decrease is reported without recoverable cell n, SE or independently paired night-dependence coefficient."
            },
            {
                "id": "giavi2020_epilobium",
                "doi": "10.1038/s41598-020-68667-y",
                "species": "Epilobium angustifolium",
                "outcome": "logit coefficient for selected capsules that produce seeds, conditional on capsule sampling",
                "exposure": "directly illuminated versus dark meadow",
                "reported_model_beta": giavi_logit,
                "reported_model_beta_SE": giavi_se,
                "derived_odds_ratio": math.exp(giavi_logit),
                "odds_ratio_Wald_95_lower": math.exp(giavi_logit-1.96*giavi_se),
                "odds_ratio_Wald_95_upper": math.exp(giavi_logit+1.96*giavi_se),
                "ln_RR": None,
                "primary_initial_flower_endpoint": False,
                "limitations": "An odds ratio is not a risk ratio. Reproduction model unit and 4 lit / 4 control meadow clustering need independent inspection."
            },
            {
                "id": "boom2020_green",
                "doi": "10.1038/s41598-020-75471-1",
                "species": "Silene latifolia",
                "outcome": "author-defined estimated seeds per originally exposed flower (product of outcome stages)",
                "exposure": "green LED versus dark",
                "derived_green_dark_ratio": 0.5150,
                "ln_RR_approx_product": math.log(0.5150),
                "ln_RR_SE": None,
                "primary_initial_flower_endpoint": False,
                "limitations": "Estimated, not observed total mature seeds; published net-seed ANOVA p=0.63; colors share experiment."
            },
        ],
        "forbidden_inferences": [
            "Do not combine four rows above into one effect-size estimate.",
            "Do not treat positive Li and negative Lu as independent species-matched slope observations.",
            "Do not turn an OR into RR or assume logistic model SE is an LnRR SE.",
            "Do not attribute all loss to moths without effective guild data and light susceptibility.",
            "Do not treat missing SE or non-significance as zero.",
        ],
    }
    return out

def write_result(path: Path) -> dict:
    obj = publication_compatible_examples()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False)+"\n",encoding="utf-8")
    return obj

if __name__=="__main__":
    target=Path("results/alan_A_source_effect_type_gate_20261115.json")
    x=write_result(target)
    print(json.dumps({"status":"verified_aggregate_arithmetic_only","example_count":len(x["studies"]),"beta_V":x["confirmed_A_beta_V"]},ensure_ascii=False))
