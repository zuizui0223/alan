#!/usr/bin/env python3
"""Counterfactual MOIST missingness sensitivity, *not* trait imputation or new ALAN biology.

Only evaluates 49 strata containing positive known MOIST cases and references.
Unknown values are never written back as observed. Hypothetical shift magnitudes
illustrate which missing-at-random assumption could reverse the selected-group contrast.
"""
import json
from scripts.analysis import alan_wet_species_support_20261117 as base


def audit():
    source,sha=base.read_source()
    strata,_=base.get_strata(source)
    observed,n_known=base.difference(strata)
    n_case=sum(g["case_all"] for g in strata.values())
    n_ref=sum(g["ref_all"] for g in strata.values())
    missing_case=sum(g["case_all"]-len(g["case"]) for g in strata.values())
    missing_ref=sum(g["ref_all"]-len(g["ref"]) for g in strata.values())
    # Expected contrast if unscored species of each role have the observed
    # mean within that role AND exact stratum, *as a mathematical scenario*.
    sumcontrast=0.
    for g in strata.values():
        ma=sum(g["case"])/len(g["case"])
        mr=sum(g["ref"])/len(g["ref"])
        sumcontrast+=g["case_all"]*(ma-mr)
    hypothetical_same_role=sumcontrast/n_case
    alpha=missing_case/n_case
    beta=sum(g["case_all"]*(g["ref_all"]-len(g["ref"]))/g["ref_all"] for g in strata.values())/n_case
    assert alpha>0 and beta>0 and hypothetical_same_role>0
    def scenario(case_shift,reference_shift):
        # Case unknown mean = observed case mean + case_shift;
        # ref unknown mean = observed ref mean + reference_shift.
        return hypothetical_same_role+alpha*case_shift-beta*reference_shift
    case_tipping=hypothetical_same_role/alpha
    ref_tipping=hypothetical_same_role/beta
    symmetric_tipping=hypothetical_same_role/(alpha+beta)
    assert abs(scenario(-case_tipping,0))<1e-10
    assert abs(scenario(0,ref_tipping))<1e-10
    assert abs(scenario(-symmetric_tipping,symmetric_tipping))<1e-10
    return {
        "schema":"alan_wet_moisture_missingness_tipping_v1",
        "source_sha256":sha,
        "scope":"hypothetical offsets within matched strata; no missing value is actually imputed",
        "source_species_total":len(source),
        "matched_original_strata":len(strata),
        "known_case_count":n_known,
        "matched_strata_case_total":n_case,
        "matched_strata_case_unknown":missing_case,
        "matched_strata_reference_total":n_ref,
        "matched_strata_reference_unknown":missing_ref,
        "observed_complete_case_difference":observed,
        "same_role_missing_equals_observed_mean_scenario":hypothetical_same_role,
        "weight_case_missing_alpha":alpha,
        "weight_reference_missing_beta":beta,
        "contrast_formula":"D = D_equal_role + alpha*shift_missing_case - beta*shift_missing_reference",
        "tipping_case_unknown_points_below_known_case":case_tipping,
        "tipping_reference_unknown_points_above_known_reference":ref_tipping,
        "tipping_symmetric_each_side_points":symmetric_tipping,
        "scenario_both_0p5_points_away":scenario(-0.5,0.5),
        "scenario_both_1p0_points_away":scenario(-1.0,1.0),
        "excluded_unmatched_strata":"The 49-stratum conditional estimand excludes 26 candidate species outside these matched strata.",
        "conclusion":"tipping thresholds quantify model dependence; there is no evidence that missing traits actually differ by these offsets"
    }


if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,sort_keys=True,indent=2))
