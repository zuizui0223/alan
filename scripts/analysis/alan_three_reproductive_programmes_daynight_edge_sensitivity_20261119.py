#!/usr/bin/env python3
"""Original published field-source coefficients and day/night pathway sensitivity.

The ONLY independently supported biological inference is the authors'
published field or controlled-lab model result, with its original units.

Macgregor 2019's four caging percentages are across lighting treatments.
A day/night route ablation is a mathematical scenario, NOT an ALAN causal
estimate. Giavi 2020's reported Table 1(c) "relative reproduction"
binomial denominator may not match the described biological denominator:
keep the beta/SE and NEVER re-label exp(beta) as measured seed RR or odds.
Wilson 2021 natural skyglow is observational and has NO verified beta/SE.
"""
import argparse
import json
import math
from pathlib import Path

SCHEMA="alan_source_3_true_plant_reproduction_pathway_contrasts_20261119_v1"
SOURCE_REVIEW_SHA="a4c105c50f7b53cb2c58b287d065b79432b6ec7c308c35b1c4b873e78a4a359e"
LEDGER=Path("results/all_species_forest/biotic_interactions/alan_three_true_plant_reproduction_published_coefficient_ledger_20261119.json")


def phi(z):
    """Standard normal cumulative; used to check only author-rounded Wald z."""
    return (1+math.erf(z/math.sqrt(2)))/2


def report(data):
    if data.get("schema")!=SCHEMA or data.get("source_review_SHA256")!=SOURCE_REVIEW_SHA:
        raise ValueError("Original 3-study provenance changed")
    scope=data.get("model_scope",{})
    if not all(scope.get(k) is True for k in (
        "reproduction_and_prey_capture_not_common_effect_size",
        "three_publications_not_three_matching_seed_effects",
        "conditional_caging_compensation_is_noncausal_sensitivity",
        "no_cross_study_random_effect_meta_fitted",
        "treatment_or_caging_group_not_equal_independent_light_installation"
    )):
        raise ValueError("Incompatible ecological stages or source independence gate changed")
    studies={s["programme"]:s for s in data["primary_studies"]}
    if set(studies)!={"Macgregor2019","Giavi2020","Wilson2021"}:
        raise ValueError("Unexpected original study identity")
    m=studies["Macgregor2019"]; g=studies["Giavi2020"]; w=studies["Wilson2021"]
    if m["doi"]!="10.1002/ecs2.2550" or g["doi"]!="10.1038/s41598-020-68667-y" or w["doi"]!="10.1093/icb/icab010":
        raise ValueError("Original DOI identity changed")
    if m["design"]["independent_field_plots"]!=5 or m["design"]["repeated_field_runs"]!=4:
        raise ValueError("Macgregor physical installation grain drift")
    if m["design"]["observed_flowers"]!=417 or m["design"]["individual_female_plants"]!=169:
        raise ValueError("Original source Macgregor flower/plant counts changed")
    if (m["design"]["pollinated_flowers"]!=194 or
        m["design"]["damaged_successfully_pollinated_capsules"]!=46 or
        m["design"]["undamaged_capsules_for_seed_count_or_mass"]!=148):
        raise ValueError("Original source damaged/seed capsule coverage changed")
    p=m["published_day_night_flower_access_success"]
    if not all(0<=p[x]<=1 for x in ("open_both","day_only","night_only")):
        raise ValueError("Invalid author published pollination success proportion")
    if abs(p["open_both"]-.72)>1e-9 or abs(p["day_only"]-.376)>1e-9 or abs(p["night_only"]-.528)>1e-9:
        raise ValueError("Original author reported caging percentages changed")
    if m["design"]["source_original_numerical_treatment_CSV_available"]:
        raise ValueError("Do not claim missing original treatment rows are available")

    day=p["day_only"];night=p["night_only"];observed=p["open_both"]
    predicted=1-(1-day)*(1-night)
    sensitivity=[]
    for pct in (0,.25,.5,.75,1):
        hypothetical_night=night*(1-pct)
        combined=1-(1-day)*(1-hypothetical_night)
        sensitivity.append({
            "hypothetical_fraction_of_night_route_success_removed":pct,
            "predicted_combined_pollination_success_if_day_unchanged":combined,
            "difference_from_independent_channels_no_night_loss_pp":100*(combined-predicted),
            "not_real_ALAN_treatment_effect":True,
            "not_fitted_visit_to_pollen_efficiency":True,
        })

    mg={
        "study_DOI":m["doi"],
        "source_original_observed_open_access_day_night_percent":{"open":observed,"day":day,"night":night},
        "independent_day_night_route_predicted_open":predicted,
        "observed_minus_independent_route_prediction_percentage_points":100*(observed-predicted),
        "published_damaged_fraction_among_successfully_pollinated":46/194,
        "damaged_not_measured_ALAN_effect":True,
        "assumed_independence_algebra_not_inferred_from_same_flower":True,
        "marginal_caging_rates_not_a_dark_control_only_comparison":True,
        "night_loss_sensitivity_not_a_causal_ALAN_dose_response":sensitivity,
        "published_full_night_pollination_direction":m["lighting_results"]["full_night_vs_dark_pollination_success"],
        "published_part_night_pollination_direction":m["lighting_results"]["part_night_vs_dark_pollination_success"],
        "source_primary_R_code_only_numeric_raw_CSV_missing":True
    }

    if g["design"]["independent_lit_meadows"]!=4 or g["design"]["independent_dark_control_meadows"]!=4:
        raise ValueError("Original Giavi lit/dark field design drift")
    if g["design"]["potted_plants"]!=97 or g["design"]["Silene_individual_plants"]!=47 or g["design"]["Epilobium_individual_plants"]!=50:
        raise ValueError("Source original Giavi plant counts drift")
    if not g["design"]["adjacent_dark_shares_treatment_meadow"]:
        raise ValueError("Adjacent dark is never a separate randomized treatment meadow")
    if g["design"]["lab_moth_replicates"]!=27 or g["design"]["lab_distinct_females"]!=22:
        raise ValueError("Original repeated lab moth grain drift")
    original_contrasts=[]
    ids=[]
    for row in g["published_model_contrasts"]:
        beta=float(row["beta"]); se=float(row["SE"])
        if not math.isfinite(beta) or not math.isfinite(se) or se<=0:
            raise ValueError("Invalid original published model beta or SE")
        key=(row["endpoint"],row["contrast"])
        if key in ids:raise ValueError("Repeated original model contrast")
        ids.append(key)
        desc={
            "study_DOI":g["doi"],"endpoint":row["endpoint"],
            "contrast":row["contrast"],"orientation":row["orientation"],
            "published_model":row["model"],
            "published_beta":beta,"published_SE":se,
            "published_p_value":row["p_published"],
            "reconstructed_link_scale_Wald_95pct_CI":[beta-1.96*se,beta+1.96*se],
            "reported_effect_scale":row["effect_type"],
            "same_field_and_plant_cluster_as_other_model_rows":True,
            "not_an_independently_replication_or_cross_study_effect":True,
        }
        if row["effect_type"]=="odds_ratio_from_reported_binomial_logit":
            desc["source_model_conditional_attack_odds_ratio"]=math.exp(beta)
            desc["attack_odds_ratio_Wald_95pct_CI"]=[math.exp(beta-1.96*se),math.exp(beta+1.96*se)]
            desc["is_not_attack_rate_risk_ratio"]=True
        elif row["effect_type"]=="original_fitted_link_beta_ONLY_due_to_documented_binomial_denominator_ambiguity":
            desc["do_not_exponentiate_as_fruit_or_seed_probability_ratio"]=True
            desc["source_reproductive_count_denominator_source_ambiguity"]=True
        else:
            desc["exponentiated_beta_not_claimed_fitness_OR_or_RR"]=True
        # Author reported p is source model; do not silently replace by a
        # calculated normal approximation from beta/SE which can differ.
        original_contrasts.append(desc)
    required={
        ("Silene_infested_fruit_fraction","adjacent dark vs independent dark control"),
        ("Silene_uninfected_capsules_per_exposed_flower","independent dark vs adjacent dark"),
        ("Epilobium_selected_capsule_at_least_four_mature_seeds","illuminated vs dark independent control"),
    }
    if not required.issubset(set(ids)):
        raise ValueError("Critical original signed table contrasts omitted")
    # NEW purely conditional seed-predation compensation threshold, using
    # ONLY the source-reported Giavi Table 1b OR.  Let p be the chance an
    # exposed flower forms a capsule, h the fraction of formed capsules
    # infected. W=p*(1-h) is the proportion of exposed flowers producing
    # at least one UNINFECTED capsule (not a seed count).
    # For h_control=h and attack odds ratio OR, h_adj=OR*h/(1-h+OR*h).
    # No net loss requires p_adj/p_control >= (1-h)/(1-h_adj)
    #   = 1+(OR-1)*h.
    # Baseline h here is HYPOTHETICAL, not read from raw original plots.
    odds_rows=[a for a in original_contrasts
               if a["endpoint"]=="Silene_infested_fruit_fraction"
               and a["orientation"]=="adjacent/dark"]
    if len(odds_rows)!=1:
        raise ValueError("Original adjacent-vs-dark attack odds must be singular")
    observed_odds=odds_rows[0]["source_model_conditional_attack_odds_ratio"]
    odds_ci=odds_rows[0]["attack_odds_ratio_Wald_95pct_CI"]
    conditional_threshold=[]
    for h in (0.05,0.10,0.20,0.30):
        infected_adj=observed_odds*h/(1-h+observed_odds*h)
        threshold=1+(observed_odds-1)*h
        conditional_threshold.append({
            "hypothetical_original_dark_fruit_infection_probability":h,
            "model_constrained_adjacent_fruit_infection_probability":infected_adj,
            "minimum_adjacent_vs_dark_capsule_formation_probability_ratio_for_no_uninfected_capsule_loss":threshold,
            "conditional_ratio_from_low_attack_OR":1+(odds_ci[0]-1)*h,
            "conditional_ratio_from_high_attack_OR":1+(odds_ci[1]-1)*h,
            "maximum_original_dark_capsule_formation_probability_where_full_rescue_is_even_possible":1/threshold,
            "NOT_an_observed_ALAN_pollination_effect":True,
            "NOT_a_total_seed_set_estimate":True,
            "baseline_infection_hypothetical_not_from_Giavi_original_samples":True,
        })
    gi={
        "study_DOI":g["doi"],
        "physical_lit_meadows":4,
        "physical_dark_control_meadows":4,
        "potted_plants_not_light_installations":97,
        "same_original_field_and_lab_experiments_not_causally_joint_fitted":True,
        "original_published_effects_separately_by_endpoint":original_contrasts,
        "cannot_pretend_published_reproduction_link_is_actual_seed_loss_ratio":True,
        "no_raw_plot_level_field_reanalysis":True,
        "conditional_compensatory_pollination_thresholds":conditional_threshold,
        "odds_ratio_to_risk_model_not_original_randomized_pollen_success":True,
    }
    wi=[]
    for row in w["published_mechanism_contrasts"]:
        if row["reported_effect_beta_SE_available"]:
            raise ValueError("Wilson verified effect SE source gate not passed")
        wi.append({
            "endpoint":row["endpoint"],
            "exposure":row["exposure"],
            "published_direction":row["direction"],
            "source_light_assignment":row["randomized"],
            "source_effect_size_beta_SE_not_verified":True,
        })
    if len(wi)!=4 or w["design"]["original_independent_light_towers_n"] is not None:
        raise ValueError("Wilson source treatment clusters must remain unverified")
    return {
        "schema":"alan_20261119_three_true_repro_published_contrast_nonmeta_sensitivity_v1",
        "original_review_source_SHA256":SOURCE_REVIEW_SHA,
        "primary_published_programmes_n":len(studies),
        "primary_direct_plant_reproduction_programmes_n":3,
        "macgregor2019_published_day_night_functional_pathways":mg,
        "giavi2020_original_field_and_lab_model_table_contrasts":gi,
        "wilson2021_distinguish_ambient_skyglow_and_direct_tower":wi,
        "primary_reproductive_effect_random_effects_meta_fitted":False,
        "three_primary_effects_with_same_endpoint_SE_and_verified_treatment_units":0,
        "cross_programme_CI_not_estimated":True,
        "unmeasured_alternative_pollination_or_selfing_cannot_be_filled_in":True,
        "spatial_edge_effect_not_equivalent_to_increased_lit_center":True,
        "original_source_primary_metadata_and_published_model_only_not_original_experiment_row_refit":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ledger",type=Path,default=LEDGER)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    result=report(json.loads(a.ledger.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    mg=result["macgregor2019_published_day_night_functional_pathways"]
    gg=result["giavi2020_original_field_and_lab_model_table_contrasts"]
    print("ALAN_MACGREGOR2019_DAY_NIGHT_INDEPENDENT_ROUTE_PRED",mg["independent_day_night_route_predicted_open"])
    print("ALAN_MACGREGOR2019_DAY_NIGHT_OBSERVED_PRED_DIFF_PP",mg["observed_minus_independent_route_prediction_percentage_points"])
    print("ALAN_GIAVI2020_ORIGINAL_TABLE_EFFECTS",[
        {"endpoint":x["endpoint"],"contrast":x["contrast"],"beta":x["published_beta"],
         "OR":x.get("source_model_conditional_attack_odds_ratio")} for x in gg["original_published_effects_separately_by_endpoint"]
    ])
    print("ALAN_3_PAPER_REPRODUCTIVE_META_GATES",result["three_primary_effects_with_same_endpoint_SE_and_verified_treatment_units"])


if __name__=="__main__":
    main()
