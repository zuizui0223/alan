#!/usr/bin/env python3
"""Reconstruct paired plot-difference SEs from Firebaugh & Haynes 2020.

Author DOI 10.1002/ecs2.3227, Results: 4 randomized lit/dark plot pairs
surveyed by sweep net on one summer date in each of 2015 and 2016.
Original lit-minus-dark density t(3) statistics have been rounded.
The *same four spatial plot pairs* were reused across the two years.

No conversion to lnRR, no claim of a significant year interaction,
no extra study k, and no substitution for predator/herbivore contrast.
"""
import argparse
import json
import math
from pathlib import Path

DOI="10.1002/ecs2.3227"
T975_DF3=3.182446305284263

ORIGINAL_AUTHOR_REPORT={
    2015:{"mean_lit_minus_dark_per_m2":6.04,"t3":2.89,"author_p_rounded":0.06,
          "author_text_lit_over_dark_ratio_approx":1.61},
    2016:{"mean_lit_minus_dark_per_m2":-0.15,"t3":-0.07,"author_p_rounded":0.95,
          "author_text_lit_over_dark_ratio_approx":None},
}
TROPHIC_INTERACTION={
    "author_term":"ALAN x Trophic group",
    "source_model":"Box-Cox transformed arthropod density (lambda=-0.16)",
    "model_coefficient":0.26,
    "model_SE":0.46,
    "model_df":24,
    "author_p_rounded":0.58,
    "not_direct_predation_rate":True,
    "not_lnRR":True,
}

def compute():
    years={}
    for year,r in ORIGINAL_AUTHOR_REPORT.items():
        d,t=r["mean_lit_minus_dark_per_m2"],r["t3"]
        if d*t<=0 or t==0:raise ValueError("STOP: source contrast sign or t discrepancy")
        se=abs(d/t)
        years[str(year)]={
            "source_reported_mean_lit_minus_dark_per_m2":d,
            "source_reported_t3":t,
            "author_reported_p_rounded":r["author_p_rounded"],
            "implied_SE_difference_from_rounded_t3":se,
            "approx_paired_t3_95_difference":[d-T975_DF3*se,d+T975_DF3*se],
            "source_reported_lit_to_dark_relative_description":r["author_text_lit_over_dark_ratio_approx"],
            "independent_spatial_pair_n":4,
            "same_four_pairs_as_other_year":True,
            "not_author_raw_pair_data":True,
        }
    d2=years["2016"]["source_reported_mean_lit_minus_dark_per_m2"]
    d1=years["2015"]["source_reported_mean_lit_minus_dark_per_m2"]
    s2=years["2016"]["implied_SE_difference_from_rounded_t3"]
    s1=years["2015"]["implied_SE_difference_from_rounded_t3"]
    change=d2-d1
    bound=[abs(s2-s1),s1+s2]
    # The two years reuse the same plots; covariance is not reported.
    # Do not use this extreme algebraic range as an estimated CI or p value.
    return {
        "doi":DOI,
        "source":"Original peer-reviewed article, Results text and Table 2",
        "physical_programme":"Blandy_2015_2016_warm_season_grassland_Firebaugh",
        "n_true_independent_randomized_plot_pairs":4,
        "n_annual_observations_per_arm":4,
        "two_annual_outcomes_are_dependently_measured":True,
        "measurement":"Above-ground sweep-net arthropod density, individuals/m^2",
        "design":"8 plots in four pairs, one independently randomly lit plot per pair, repeated in 2015 and 2016",
        "source_rounded_paired_differences":years,
        "difference_2016_minus_2015_in_pair_mean_differences":change,
        "unknown_year_to_year_covariance_SE_bounds":[min(bound),max(bound)],
        "cannot_infer_year_by_light_interaction_p_without_joint_covariance":True,
        "author_trophic_group_interaction":TROPHIC_INTERACTION,
        "author_reported_total_collected":9925,
        "do_not_use_arthropod_individuals_as_independent_ALAN_units":True,
        "quantitative_meta_eligible_today":False,
        "why_not_in_stage_meta":"Paired mean density difference not lnRR; original 4-pair lit/dark values or SE of log response, year covariance, and trophic-group coeffects required for multi-study comparable estimand.",
        "independent_meta_programmes_added":0,
        "status":"AUTHOR_T_ROUNDED_PAIRED_DIFF_SE_RECONSTRUCTED_NOT_HARMONIZED_META"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("results/alan_firebaugh2020_paired_t_audit_20261124.json"))
    x=p.parse_args()
    result=compute()
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("FIREBAUGH2020_ROUNDED_AUTHOR_T",json.dumps({
        "2015":result["source_rounded_paired_differences"]["2015"],
        "2016":result["source_rounded_paired_differences"]["2016"],
        "year_change_point":result["difference_2016_minus_2015_in_pair_mean_differences"],
        "year_change_SE_bounds":result["unknown_year_to_year_covariance_SE_bounds"],
        "meta_k_increment":0},ensure_ascii=False))

if __name__=="__main__":
    main()
