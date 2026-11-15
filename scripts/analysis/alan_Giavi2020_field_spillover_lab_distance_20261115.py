#!/usr/bin/env python3
"""New original-source Table2 Hadena lit/dark distance gradients tied to
previously established Giavi2020 Table1 fruit-infestation dark spillover.

Source original JATS DOI 10.1038/s41598-020-68667-y, PMC7368033.
Table1 field predation OR existed in alan already (2026-11-14).
Table2 light/distance moth visits and egg laying original model coefficients
have NOT yet been recorded in the same source-verified pathway ledger.

One paper with 8 field meadows (4 new-light and 4 controls) and
27 lab moth release sessions, 22 individual females. Lit and adjacent
dark occupy the SAME four treated meadows. Exact 27 sessions are not 216
independently randomized moths. No actual stigma pollen or raw plants.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from scripts.analysis.alan_crossstudy_scope_20261114 import acquire_giavi_tables,require
from scripts.analysis.alan_giavi2020_reported_contrasts_20261114 import sections

DOI="10.1038/s41598-020-68667-y"
EXPECTED={
 "field_seedpred_dark_vs_adj":(-1.90,.576),
 "field_seedpred_lit_vs_adj":(-1.14,.363),
 "field_unattacked_repro_dark_vs_adj":(1.04,.499),
 "field_unattacked_repro_lit_vs_adj":(.69,.371),
 "field_epil_lit_seed_presence":(-1.37,.726),
 "lab_visit_dark_dist_lamp":(-.08,.075),
 "lab_visit_lit_dist_lamp":(.37,.104),
 "lab_visit_lit_intercept":(-1.22,.522),
 "lab_eggs_dark_dist_lamp":(.06,.094),
 "lab_eggs_lit_dist_lamp":(.24,.123),
 "lab_eggs_lit_intercept":(-1.17,.622),
 "lab_total_eggs_lit":(.15,.239)}

def checked(row,tag):
    beta,se=EXPECTED[tag]
    require(math.isclose(row["beta"],beta,abs_tol=1e-8) and
            math.isclose(row["se"],se,abs_tol=1e-8),
            "Giavi2020 authentic original Table coefficient drift "+tag)
    return row

def transform(beta,se,link):
    assert link in ("logit","log")
    z=1.959963984540054
    return {"author_model_beta":beta,"author_SE":se,"model_link":link,
       "exponentiated_coefficient":math.exp(beta),
       "Wald95_from_rounded_published_SE":[math.exp(beta-z*se),math.exp(beta+z*se)],
       "exponentiation_is_odds_ratio_for_logit_rate_multiplier_for_log":True,
       "not_original_plantlevel_fitted_regression":True}

def analyze(original):
    require(original["article_doi"]==DOI,"Wrong source DOI")
    tables={x["label"]:x for x in original["original_published_model_tables"]}
    require(set(tables)>={"Table 1","Table 2"},"Need both primary JATS tables")
    a=sections(tables["Table 1"]["rows"])
    b=sections(tables["Table 2"]["rows"])
    def extract(group,label,tag,link):
        row=checked(group[label],tag)
        return {**transform(row["beta"],row["se"],link),
                 "original_published_p":row["reported_p"],"original_source_label":label}
    field={
      "adjacent_dark_vs_far_dark_infected_fruit_odds":{
        **transform(-checked(a["(b)"]["Dark"],"field_seedpred_dark_vs_adj")["beta"],
                    a["(b)"]["Dark"]["se"],"logit"),
        "original_author_field_p":a["(b)"]["Dark"]["reported_p"],
        "contrast":"Adjacent dark / far dark",
        "previously_recomputed_in_ALAN_Nov14":True},
      "adjacent_dark_vs_direct_lit_infected_fruit_odds":{
        **transform(-checked(a["(b)"]["Illuminated"],"field_seedpred_lit_vs_adj")["beta"],
                    a["(b)"]["Illuminated"]["se"],"logit"),
        "original_author_field_p":a["(b)"]["Illuminated"]["reported_p"],
        "contrast":"Adjacent dark / direct light",
        "previously_recomputed_in_ALAN_Nov14":True},
      "unattacked_fruit_reproduction_dark_vs_adjacent_source_link":{
        "beta":checked(a["(c)"]["Dark"],"field_unattacked_repro_dark_vs_adj")["beta"],
        "se":a["(c)"]["Dark"]["se"],
        "original_p":a["(c)"]["Dark"]["reported_p"],
        "response":"proportion of exposed flowers becoming unattacked fruit capsules",
        "source_bivariate_binomial_matrix_denominator_ambiguous":True,
        "do_not_exponentiate_as_confirmed_true_plant_success_odds":True},
      "unattacked_fruit_reproduction_lit_vs_adjacent_source_link":{
        "beta":checked(a["(c)"]["Illuminated"],"field_unattacked_repro_lit_vs_adj")["beta"],
        "se":a["(c)"]["Illuminated"]["se"],
        "original_p":a["(c)"]["Illuminated"]["reported_p"],
        "source_bivariate_binomial_matrix_denominator_ambiguous":True},
      "Epilobium_lit_vs_dark_selected_capsule_seed_presence_odds":
        extract(a["(a)"],"Illuminated","field_epil_lit_seed_presence","logit")}
    visit={
       "dark_distance_per_1m":extract(b["(a)"],"Dark: distance from lamp",
                                     "lab_visit_dark_dist_lamp","log"),
       "lit_distance_per_1m":extract(b["(a)"],"Illuminated: distance from lamp",
                                    "lab_visit_lit_dist_lamp","log"),
       "lit_conditional_at_model_zero_distances":extract(
                       b["(a)"],"Illuminated","lab_visit_lit_intercept","log")}
    egg={
       "dark_distance_per_1m":extract(b["(b)"],"Dark: distance from lamp",
                                     "lab_eggs_dark_dist_lamp","log"),
       "lit_distance_per_1m":extract(b["(b)"],"Illuminated: Distance from lamp",
                                    "lab_eggs_lit_dist_lamp","log"),
       "lit_conditional_at_model_zero_distances":extract(
                       b["(b)"],"Illuminated","lab_eggs_lit_intercept","log")}
    for d in (visit,egg):
        for name in ("dark_distance_per_1m","lit_distance_per_1m"):
            d[name]["source_two_meter_distance_multiplier"]=math.exp(
                 2*d[name]["author_model_beta"])
            d[name]["two_meter_multiplier_does_NOT_measure_individual_visit_probability"]=True
    total=extract(b["(c)"],"Illuminated","lab_total_eggs_lit","log")
    require(visit["lit_distance_per_1m"]["author_model_beta"]>0 and
            egg["lit_distance_per_1m"]["author_model_beta"]>0,
            "Lab expected spatial direction absent")
    require(visit["dark_distance_per_1m"]["author_model_beta"]<0,
            "Lab dark original location slope changed")
    require(field["adjacent_dark_vs_far_dark_infected_fruit_odds"]["exponentiated_coefficient"]>6,
            "Previously documented original field spatial predation absent")
    return {
     "schema":"ALAN_Giavi2020_original_Table1_field_spillover_and_Table2_lab_moth_distance_v1",
     "original_paper_DOI":DOI,"PMCID":"PMC7368033",
     "reacquired_primary_JATS_sha256":original["article_xml_sha256"],
     "published_field_design":{"true_original_meadow_clusters":8,
       "distinct_artificially_lit_meadows":4,
       "distinct_far_dark_control_meadows":4,
       "adjacent_dark_and_direct_light_share_same_four_meadows":True,
       "total_potted_phytometer_plants":97,
       "S_latifolia_plants":47,"E_angustifolium_plants":50,
       "E_angustifolium_selected_capsules_seed_present_count_assessed":240,
       "S_latifolia_selected_uninfested_capsules_count_assessed":199,
       "mean_lux_adjacent_dark":.0021,"mean_lux_far_dark_control":.0040,
       "mean_lux_illuminated":16.1,
       "adjacent_ambient_lux_lower_than_far_dark_yet_higher_seed_predation":True},
     "published_lab_design":{"independent_moth_individuals":22,
       "true_release_sessions":27,
       "lit_sessions":16,"dark_sessions":11,
       "flowers_checked_for_eggs":216,
       "lab_corridor_structures":3,
       "flower_distance_to_lamp_meters":[0,2,4,6],
       "original_response_models_include_release_distance_and_moth_random_effect":True,
       "reuse_of_some_moths_between_sessions":True},
     "original_field_stage_coefficients":field,
     "original_lab_Table2_visitation_response":visit,
     "original_lab_Table2_oviposition_response":egg,
     "original_lab_Table2_total_egg_count_conditional_light_effect":total,
     "source_inference_limits":{
       "Table1_spillover_Beta_OR_was_ALREADY_IN_repo_20261114":True,
       "NEW_Table2_light_distance_gradient_with_SE_reexpressed":True,
       "moth_distance_slope_difference_SE_unavailable_without_model_covariance":True,
       "in_lit_corridor_moth_visits_and_oviposition_increase_away_from_lamp":True,
       "but_lab_field_sequential_link_NOT_actual_mediation_estimation":True,
       "field_and_lab_from_ONE_PUBLICATION_not_two_independent_replicates":True,
       "female_moth_visits_not_measured_effective_conspecific_pollen":True,
       "seed_predation_from_moth_LARVAE_not_moth_predator_gecko_bat":True,
       "adjacent_dark_site_light_exposure_not_local_lux_only":True,
       "relative_reproduction_original_denominator_syntax_ambiguous":True,
       "no_new_individuallevel_repro_model_fitted":True,
       "no_island_mainland_comparison_or_global_pooled_coefficient":True}}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();output=a.out
    require(not output.exists(),"Do not overwrite source artifacts")
    z=analyze(acquire_giavi_tables())
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(z,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print("ALAN_GIAVI2020_OLD_FIELD_SPILLOVER_NEW_LAB_DISTANCE",json.dumps({
      "source":z["reacquired_primary_JATS_sha256"],
      "field":z["original_field_stage_coefficients"],
      "lab_visits":z["original_lab_Table2_visitation_response"],
      "lab_eggs":z["original_lab_Table2_oviposition_response"],
      "no_mediation_claimed":True},ensure_ascii=False))
if __name__=="__main__":main()
