#!/usr/bin/env python3
"""Exploratory original Figure-4 plot-readout: ALAN phenological opportunity.

Source: Grognuz, Gisler & Knop (2026), 14 light/dark paired site-years, 16
wildflowers, 4 stage-specific mixed-effects author models. Values here are
approximate figure diamond positions, NOT authors' original regression table
or plot-level data. All stage estimates within one study are dependent.

No co-flowering overlap or plant-pollinator mismatch is measured: without
absolute dark-site phenology or visitors, only DIFFERENTIAL treatment
responses can be described. Sensitivity ±3 days is illustrative figure-reading
perturbation, not a statistical 95% interval or original author uncertainty.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
from statistics import mean, median

DEFAULT=Path("results/all_species_forest/biotic_interactions/grognuz2026_figure4_16plant_phenology_visual_approx_source_20261119.json")
STAGES=("FIRST","QUARTER","PEAK","ULTIMATE")
SCHEMA="alan_grognuz2026_original_figure4_approx_treatment_shift_day_source_register_v1"


def analyze(registry):
    if registry.get("schema")!=SCHEMA or registry.get("primary_doi")!="10.1016/j.biocon.2025.111650":
        raise ValueError("Original published Figure 4 source identity has drifted")
    d=registry["design"]
    if d["matched_light_dark_site_pairs_in_2022"]!=6 or d["matched_light_dark_site_pairs_in_2023"]!=8:
        raise ValueError("2022/2023 independent lighting pair-year designs differ from original")
    if d["matched_site_pair_years_total"]!=14 or d["within_site_inner_pheno_subplots"]!=108:
        raise ValueError("Original 14 pair-years and 108 nested plots drifted")
    if d["actual_pollinator_visitation_or_seed_fitness_in_same_experiment_measured"]:
        raise ValueError("This study DOES NOT identify realized pollinator overlap or seeds")
    figure=registry["source_figure_4"]
    if not figure["not_original_data_or_published_exact_coefficients"]:
        raise ValueError("Visual readout must never be passed as original exact coefficients")
    if not figure["visual_readout_uncertainty_unquantified"]:
        raise ValueError("Source figure-readout error is not a measured covariance")
    stageinfo=figure["source_species"]
    if len(stageinfo)!=16 or d["original_analyzed_flowering_species"]!=16:
        raise ValueError("16 original flowering species identity drift")
    if d["original_published_species_with_at_least_one_week_shift_at_any_stage"]!=12:
        raise ValueError("Author original 12/16 species condition drift")
    names=[x["species"] for x in stageinfo]
    if len(set(names))!=16:
        raise ValueError("Duplicate genus/species would pseudo-replicate source readout")
    e=int(figure["illustrative_artificial_rounding_perturbation_days"])
    if not 0<=e<=10:
        raise ValueError("A different figure transcription perturbation requires review")
    if not figure["readout_uncertainty_plus_minus_3_days_is_not_source_confidence_interval"]:
        raise ValueError("Illustrative bounds must not be mistaken for 95% source confidence interval")
    for x in stageinfo:
        if set(x)!={"species",*STAGES}:
            raise ValueError("Unrecognized author plotted phenological stage")
        if not all(isinstance(x[y],(int,float)) and math.isfinite(x[y])
                   and abs(x[y])<=50 for y in STAGES):
            raise ValueError("Invalid figure-readout approximate day shift")
    fieldstage={}
    for s in STAGES:
        values=[x[s] for x in stageinfo]
        fieldstage[s]={
            "observed_fig4_approx_min_days":min(values),
            "observed_fig4_approx_max_days":max(values),
            "approx_species_median_of_plotted_diamonds":median(values),
            "approx_species_mean_of_plotted_diamonds":mean(values),
            "source_estimates_not_independent_light_experiments":True,
            "source_statistical_CI_not_available_for_these_plot_readouts":True,
        }
    durations=[{"species":x["species"],"approx_window_duration_change_days":x["ULTIMATE"]-x["FIRST"]}
               for x in stageinfo]
    deltas=[r["approx_window_duration_change_days"] for r in durations]
    # If visual median for each of FIRST and LAST is within ±e of the
    # original published value, duration difference error could be ±2e.
    # This is an analyst-declared thought experiment, NOT proven OCR accuracy.
    strong_extend=[x["species"] for x in durations if x["approx_window_duration_change_days"]>10+2*e]
    strong_contract=[x["species"] for x in durations if x["approx_window_duration_change_days"]<-(10+2*e)]
    directional_reversal={
        "first_earlier_and_last_later_fig4_diamond_only":[
            x["species"] for x in stageinfo if x["FIRST"]<0 and x["ULTIMATE"]>0],
        "first_later_and_last_earlier_fig4_diamond_only":[
            x["species"] for x in stageinfo if x["FIRST"]>0 and x["ULTIMATE"]<0],
    }
    # 16 choose 2 =120 pairwise differences in treatment PEAK effects, not
    # actual treated-vs-dark pairwise co-flowering overlap.
    raw_diffs=[abs(a["PEAK"]-b["PEAK"]) for a,b in itertools.combinations(stageinfo,2)]
    if len(raw_diffs)!=120:
        raise ValueError("Original 16-flower cross-species comparison count changed")
    threshold=14
    central=sum(difference>=threshold for difference in raw_diffs)
    guaranteed_under_visual_perturbation=sum(difference>=threshold+2*e for difference in raw_diffs)
    possibly_over_under_visual_perturbation=sum(difference>=max(0,threshold-2*e) for difference in raw_diffs)
    # Preserved authors' exact published outcome of daily petal movement,
    # scientifically distinct from 16 species seasonal chart approximations.
    petal=d["published_petals_shift_example"]
    if (petal["day_flowering_Malva_moschata_morning_opening_min"]!=-147
        or petal["night_flowering_Silene_noctiflora_morning_begin_closing_min"]!=80):
        raise ValueError("Published exact diel petal estimates changed")
    return {
        "schema":"alan_grognuz2026_original_published_16plant_phenological_opportunity_response_approx_v1",
        "original_source_publication":registry["primary_publication"],
        "original_published_DOI":registry["primary_doi"],
        "evidence":{
            "original_programme_count":1,
            "independent_site_pair_years":14,
            "different_physical_pairs_across_two_years_unverified":True,
            "nested_original_inner_phenology_plots":108,
            "source_flowering_species":16,
            "author_reported_week_scale_shift_species":12,
            "pollen_or_pollinator_realized_overlap_not_measured":True,
            "original_source_Mendeley_raw_data_verified":False,
            "published_figure_plot_readout_approximate_not_original_raw_model":True,
        },
        "figure4_approx_by_flowering_stage":fieldstage,
        "figure4_approx_bloom_window_changes":{
            "species_median_of_LAST_minus_FIRST_days":median(deltas),
            "species_mean_of_LAST_minus_FIRST_days":mean(deltas),
            "observed_min_max_days":[min(deltas),max(deltas)],
            "number_of_species_with_plot_shift_duration_lengthening_approx_positive":sum(x>0 for x in deltas),
            "number_of_species_with_plot_shift_duration_shortening_approx_negative":sum(x<0 for x in deltas),
            "guaranteed_more_than_10d_extension_only_IF_visual_approx_within_declared_error_days":strong_extend,
            "guaranteed_more_than_10d_contraction_only_IF_visual_approx_within_declared_error_days":strong_contract,
            "individual_species_approx_delta_duration_days":durations,
            "source_joint_CI_of_FIRST_minus_LAST_not_available":True,
            "not_a_fitted_duration_outcome":True,
            "not_actual_coflowering_length_or_fitness":True,
        },
        "figure4_approx_interspecies_PEAK_response_difference":{
            "distinct_species_pair_differences":120,
            "threshold_absolute_differential_shift_days":threshold,
            "approx_pairs_differential_shift_at_least_threshold":central,
            "minimum_guaranteed_only_IF_fig4_readout_within_illustrative_perturbation":guaranteed_under_visual_perturbation,
            "maximum_possible_only_IF_fig4_readout_within_illustrative_perturbation":possibly_over_under_visual_perturbation,
            "plotted_approx_effect_range_days":[min(x["PEAK"] for x in stageinfo),max(x["PEAK"] for x in stageinfo)],
            "NOT_true_overlap_or_relative_flowering_date":True,
            "120_pairs_not_independent_experiments":True,
        },
        "figure4_approx_directional_FIRST_ULTIMATE_reversals":directional_reversal,
        "author_original_diel_petals":{
            "Malva_moschata_day_flowering_morning_opening_min":petal["day_flowering_Malva_moschata_morning_opening_min"],
            "Malva_original_95pct_CI":petal["Malva_author_reported_CI_minutes"],
            "Silene_noctiflora_morning_start_closing_min":petal["night_flowering_Silene_noctiflora_morning_begin_closing_min"],
            "Silene_original_95pct_CI":petal["Silene_author_reported_CI_minutes"],
            "only_Malva_CI_excludes_zero":True,
            "diel_opening_not_pollen_deposition_or_seed_set":True,
        },
        "gates":{
            "16_species_same_14_site_pair_year_source_not_16_independent_experiments":True,
            "4_phenology_stage_effects_per_species_correlated":True,
            "±3d_readout_perturbation_illustrative_not_original_confidence_interval":True,
            "no_original_plot_level_SE_or_species_covariance":True,
            "not_plant_pollinator_realized_mismatch":True,
            "no_p_value_or_global_ALAN_mechanistic_moderator_fitted":True,
            "no_global_phenology_source_meta_effect_fitted":True,
        }
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=DEFAULT)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=analyze(json.loads(a.source.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_GROGNUZ_2026_FIGURE4_SOURCE_STAGE_SUMMARIES",z["figure4_approx_by_flowering_stage"])
    print("ALAN_GROGNUZ_2026_FIGURE4_POTENTIAL_PHENO_REARRANGEMENT",z["figure4_approx_interspecies_PEAK_response_difference"])
    print("ALAN_GROGNUZ_2026_FIGURE4_WINDOW_CHANGES",z["figure4_approx_bloom_window_changes"])

if __name__=="__main__":
    main()
