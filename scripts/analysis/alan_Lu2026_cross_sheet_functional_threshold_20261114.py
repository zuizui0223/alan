#!/usr/bin/env python3
"""Source-original Lu 2026 two-worksheet reconciliation and conditional visit-service sensitivity.

This is a diagnostic of two DIFFERENT original workbook sheets, not a pollen
replacement model. All results are aggregated; original workbook is SHA-gated
and no plant/block-level raw observations leave the workflow.
"""
from __future__ import annotations
import argparse
import io
import json
from pathlib import Path
from openpyxl import load_workbook
from scripts.analysis.alan_original_Lu2026_daynight_guild_syndrome_source_20261114 import (
    XLSX_SHA, download, original_design_model, require, sheet, species_descriptive,
)

SCENARIO_EFFICIENCY = (0.5, 1.0, 2.0, 3.0, 5.0, 10.0)

def read_counts(controls, source):
    result = {}
    for short, field in (("moth","Lepidoptera"),("fly","Diptera"),("total","Total")):
        c = source[field]
        result[short] = {
            "dark": float(c["original_dark_total_visit_events"]),
            "lit": float(c["original_lit_total_visit_events"]),
        }
    return result

def plant_counts(source):
    vals = source["treatment_visit_counts_UNADJUSTED_FOR_TRUE_FLOWER_AVAILABILITY"]
    return {
        "plants": int(source["plant_taxa_with_both_lit_unlit_positive_visit_records"]),
        "moth": {"dark": float(vals["moth"]["unlit_source_positive_visit_counts"]),
                 "lit": float(vals["moth"]["lit_source_positive_visit_counts"])},
        "fly": {"dark": float(vals["fly"]["unlit_source_positive_visit_counts"]),
                "lit": float(vals["fly"]["lit_source_positive_visit_counts"])},
        "total": {"dark": float(vals["raw_visits"]["unlit_source_positive_visit_counts"]),
                  "lit": float(vals["raw_visits"]["lit_source_positive_visit_counts"])},
    }

def weighted_two_order_ratio(counts, relative_fly_per_visit):
    """Pure scenario: 1 moth-like visit carries unit function; fly carries r units."""
    r = float(relative_fly_per_visit)
    if r < 0:
        raise ValueError("Relative per-visit function must be nonnegative")
    dark = counts["moth"]["dark"] + r * counts["fly"]["dark"]
    lit = counts["moth"]["lit"] + r * counts["fly"]["lit"]
    return lit / dark if dark > 0 else None

def cohort_sensitivity(counts):
    lost_moth = counts["moth"]["dark"] - counts["moth"]["lit"]
    gained_fly = counts["fly"]["lit"] - counts["fly"]["dark"]
    threshold = lost_moth / gained_fly if lost_moth > 0 and gained_fly > 0 else None
    return {
        "source_counts": counts,
        "moth_visit_loss_if_positive": lost_moth,
        "fly_visit_gain_if_positive": gained_fly,
        "fly_gain_divided_by_moth_loss_visits_only":
            gained_fly / lost_moth if lost_moth > 0 else None,
        "relative_fly_to_lepidoptera_per_visit_effectiveness_break_even_two_orders_only":
            threshold,
        "hypothetical_two_order_lit_to_dark_ratio_by_relative_fly_efficiency": [
            {"fly_relative_to_lepidoptera_efficiency":r,
             "hypothetical_service_ratio_two_orders_ONLY":weighted_two_order_ratio(counts,r)}
            for r in SCENARIO_EFFICIENCY
        ],
        "effective_pollen_or_fruit_measured":False,
        "other_insect_orders_and_true_open_flower_denominator_NOT_INCLUDED":True,
    }

def analyze(surveys, plant):
    period = surveys["Nocturnal"]
    g = period["experimental_survey_pair_group_effects"]
    source = {
        "moth":{"dark":float(g["Lepidoptera"]["original_dark_total_visit_events"]),
                "lit":float(g["Lepidoptera"]["original_lit_total_visit_events"])},
        "fly":{"dark":float(g["Diptera"]["original_dark_total_visit_events"]),
               "lit":float(g["Diptera"]["original_lit_total_visit_events"])},
        "total":{"dark":float(g["Total"]["original_dark_total_visit_events"]),
                 "lit":float(g["Total"]["original_lit_total_visit_events"])}
    }
    docs = plant["by_period"]["Nocturnal"]
    groups = docs["realized_control_moth_visitation_majority_groups"]
    majority = plant_counts(groups["source_unlit_moth_fraction_above_half"])
    other = plant_counts(groups["source_unlit_moth_fraction_not_over_half"])
    sums = {
        key: {arm:majority[key][arm]+other[key][arm] for arm in ("dark","lit")}
        for key in ("moth","fly","total")
    }
    # If these were a literal, same-unit subset of the 400 survey-level
    # totals, no observed focal insect-order tally could EXCEED an all-survey tally.
    violations = [
        {"guild":key,"treatment":arm,"focal_total":sums[key][arm],
         "all_survey_total":source[key][arm],
         "focal_minus_all_survey":sums[key][arm]-source[key][arm]}
        for key in ("moth","fly","total") for arm in ("dark","lit")
        if sums[key][arm] > source[key][arm] + 1e-9
    ]
    transfer = docs["cross_block_control_only_visitation_syndrome_transfer"]
    result = {
        "schema":"ALAN_original_Lu2026_two_worksheet_reconciliation_functional_sensitivity_20261114_v1",
        "source":"Lu et al. Ecology and Evolution 2026 DOI 10.1002/ece3.74368",
        "source_original_workbook_SHA256":XLSX_SHA,
        "source_original_survey_rows":400,
        "source_original_focal_plant_rows":1038,
        "independent_experimental_blocks":5,
        "analysis_status":"ORIGINAL_SOURCE_SHEETS_NOT_STRICTLY_NESTED_POLLEN_EFFECTIVENESS_NOT_MEASURED",
        "cross_sheet":{
            "all_400_survey_counts":source,
            "focal_11_plant_sums":sums,
            "focal_11_plant_sum_exceeds_survey_total":bool(violations),
            "strict_same_unit_focal_subset_contradiction":violations,
            "interpretation":"Original sheets cannot be pooled as one common count partition; different observation counting or worksheet definitions require author-method reconciliation. This is NOT proof the raw paper is wrong.",
            "never_subtract_focal_counts_from_survey_as_unobserved_species":True
        },
        "within_sheet_visit_only_sensitivity":{
            "400_survey_nocturnal":cohort_sensitivity(source),
            "1038_plant_rows_focal_5_moth_majority":cohort_sensitivity(majority),
            "1038_plant_rows_focal_6_other":cohort_sensitivity(other),
            "scenario_def":"One recorded Lepidoptera visit equals 1 hypothetical function unit, Diptera visit equals r units, other guilds omitted; r* is algebraic break-even, not measured pollen efficacy."
        },
        "generalization_gate":{
            "genuinely_heldout_field_blocks_eligible":int(transfer["n_supported_heldout_blocks"]),
            "eligible_blocks_majority_stronger_visitation_decline":int(
                transfer["heldout_blocks_where_moth_majority_stronger_relative_visitor_loss"]),
            "pooled_5_vs_6_classification_transfers_reliably":False,
            "effective_moth_dependent_floral_syndrome_established":False
        },
        "limits":[
            "Lepidoptera includes moths and butterflies, not certified effective nocturnal moth pollinators",
            "No stigma pollen, cross-species per-visit efficacy, flower availability, fruit or seed in this original experiment",
            "One Tibet experiment is not a global multi-study test or a verified island-mainland pollinator shift",
            "Full 400-survey total abundance differs from sum of insect order counts in 23 records",
            "Focal 11-plant original source rows do not partition the 400-survey guild counts",
            "Plant strata are unbalanced observational subgroups and independent-block transfer fails in half eligible blocks"
        ]
    }
    return result

def main(outpath):
    outpath = Path(outpath)
    require(not outpath.exists(),"Refuse overwrite")
    raw,_ = download()
    book = load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    try:
        survey = original_design_model(sheet(book,"visitation"))
        species = species_descriptive(sheet(book,"species visitation"))
    finally:
        book.close()
    result = analyze(survey,species)
    x=result["cross_sheet"]
    require(x["focal_11_plant_sum_exceeds_survey_total"],
            "Source worksheets unexpectedly nest; re-evaluate original protocol")
    require(result["generalization_gate"]["genuinely_heldout_field_blocks_eligible"]==4,
            "Source valid-block transfer gate changed")
    require(result["generalization_gate"]["eligible_blocks_majority_stronger_visitation_decline"]==2,
            "Source valid-block contrast changed")
    outpath.parent.mkdir(parents=True,exist_ok=True)
    outpath.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_LU2026_CROSS_SHEET_AND_SENSITIVITY_AUDIT",
          json.dumps(result,ensure_ascii=False),flush=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    a = p.parse_args()
    main(a.out)
