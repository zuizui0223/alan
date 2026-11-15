#!/usr/bin/env python3
"""Two independent original ALAN plant reproductive experiments, source-locked
stage identity comparison without pretending their endpoints are exchangeable.

Li et al 2026 (bioRxiv preprint): plant seed/initial flower buds =
  fruits/initial buds * seeds/fruit. 3 actual source fields. Original
  author 79 rows SHA-pinned through Github.
Boom et al 2020 (Scientific Reports): ESTIMATED seeds/open flower =
  fruit set * (1 - moth-larval fruit infestation | fruit)
    * (ovules fertilized | uninfested fruit) * 404 ovules.
  Four true original illuminated field sites; source aggregate Table 1
  via genuine XML DOI and exact printed numbers. Not independent seed data.
No observation-level causal pathway mediation. No pooled two-study meta effect.
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.request
from pathlib import Path
import xml.etree.ElementTree as ET
from scripts.analysis.alan_Li2026_real_geometry_seed_hayfield_source_20261114 import (
    source_bytes,parse,GIT_BLOB_SHA1
)
from scripts.analysis.alan_Li2026_real_two_stage_seed_fruit_20261115 import analyze as li_analyze
from scripts.analysis.alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114 import (
    XML_URL,ARTICLE_DOI,FIELD_SITES,EXPECTED,SOURCE_NAMES,
    get_table,decomposition,require
)

def log_factors(values):
    if not values or any(v<=0 for v in values.values()):
        raise ValueError("positive multiplicative stages required for log identity")
    v={key:math.log(float(value)) for key,value in values.items()}
    return {"stage_ratio":values,"log_stage_ratio":v,"sum_log_stage":sum(v.values()),
            "product_stage_ratio":math.prod(values.values())}

def boom_original():
    req=urllib.request.Request(XML_URL,headers={
      "User-Agent":"ALAN-true-two-independent-reproductive-study-stage-source/1.0",
      "Accept":"application/xml,text/xml,*/*"})
    with urllib.request.urlopen(req,timeout=45) as resp:
        blob=resp.read(5_000_001)
    require(len(blob)<5_000_001,"Primary publisher XML too large to audit")
    root=ET.fromstring(blob)
    table=get_table(root)
    dec=decomposition(table)
    source={
      "DOI":ARTICLE_DOI,"PMCID":"PMC7591485",
      "source_JATS_url":XML_URL,
      "source_fulltext_sha256":hashlib.sha256(blob).hexdigest(),
      "true_field_sites_original":FIELD_SITES,
      "published_original_Table_1_values_verbatim_source_audited":True,
      "published_original_estimated_seeds_per_flower_omnibus_ANOVA_p":0.63,
      "fruit_set_omnibus_test_p":0.17,
      "fruit_infestation_test_p_less_than":0.001,
      "ovule_fertilization_test_p_less_than":0.01,
      "seed_output_is_derived_from_source_means_not_individual_seed_count":True}
    contrasts={}
    for name in ("White","Green","Red"):
        c=dec[name]
        part=log_factors(c["algebraic_source_product_ratio_factors"])
        assert math.isclose(part["product_stage_ratio"],
                            c["ratio_to_dark_recomputed"],rel_tol=1e-10)
        contrasts[name]={
          "light_spectrum":name,
          "source_dark_reference":"Dark",
          "source_light_arm":c,
          "stage_accounting":part,
          "reported_aggregate_seed_output_test_not_significant_p063":True,
          "not_a_source_independent_seed_ovule_replicate":True}
    return {"source":source,"source_treatment_contrasts":contrasts}

def compare():
    li=li_analyze(parse(source_bytes()))
    boom=boom_original()
    li_contrast=li["linear_vs_dark_observed_composite"]
    log_bud=log_factors({
      "fruit_formation_per_initial_flower_bud":
          li_contrast["observed_fruits_per_initial_bud_ratio"],
      "seed_count_per_recorded_fruit":
          li_contrast["observed_seeds_per_fruit_ratio"]})
    assert math.isclose(log_bud["product_stage_ratio"],
                        li_contrast["observed_seeds_per_initial_bud_ratio"],
                        rel_tol=1e-12)
    li_valid={
      "paper_DOI":"10.64898/2026.07.10.737511",
      "preprint_not_confirmed_peerreviewed":True,
      "author_original_sha1_git_blob":GIT_BLOB_SHA1,
      "true_source_field_blocks":3,
      "actual_source_plants":79,"known_seed_outcome_plants":78,
      "observed_array_vs_dark_seed_per_initial_bud_stage_accounting":log_bud,
      "observed_all_three_species_Fruit_bud_ratio_positive":True,
      "within_taxon_seeds_per_fruit_L_over_dark":{
        s:q["linear_vs_dark_source_stage_ratios"]["observed_seeds_per_fruit_ratio"]
        for s,q in li["within_taxon_actual_reproductive_stage_ratios"].items()},
      "source_pooled_seed_per_fruit_L_over_dark_positive_but_without_ZC_negative":True,
      "not_causal_intermediate_or_moth_pollen_measured":True}
    b=boom["source_treatment_contrasts"]
    assert 0.50<b["Green"]["stage_accounting"]["product_stage_ratio"]<0.53
    assert li_valid["observed_array_vs_dark_seed_per_initial_bud_stage_accounting"]["product_stage_ratio"]>4
    return {
      "schema":"ALAN_P3_two_authentic_experiments_seed_fruit_seedpredation_pathway_stage_identity_v1",
      "source_independent_publications":2,
      "source_independent_field_sites_per_paper":{"Li2026":3,"Boom2020":4},
      "true_light_treatments_comparable_as_same_exposure":False,
      "true_reproductive_endpoints_comparable_as_same_measure":False,
      "source_Li2026_reproductive_stages":li_valid,
      "source_Boom2020_reproductive_stages":boom,
      "cross_study_stage_logic":{
        "Li2026_seeds_per_initial_bud":"fruit/initial bud × seeds/fruit",
        "Boom2020_estimated_seeds_per_flower":"fruit set × (1-larval infestation | fruit) × fertilized ovules per uninfested fruit × 404 original ovules",
        "Li_higher_seed_per_bud_from_more_fruit_opportunity":True,
        "Boom_color_changes_in_seedpredation_and_ovule_success_even_if_fruit_set_stable":True,
        "do_not_equate_Li_first_bud_to_Boom_source_mature_flower":True,
        "do_not_ascribe_single_statistically_significant_net_seed_compensation":True,
        "no_meta_analytic_summary_effect_estimated":True,
        "no_ADAPTIVE_shift_or_island_moderator_estimated":True,
        "no_same_species_plant_breeding_system_x_light_factorial_test":True,
        "Boom_seed_predator_Hadena_eggs_larvae_not_same_as_adult_moth_killing_gecko":True,
        "reported_Boom_seeds_per_flower_p063_not_proof_of_equivalence":True}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True)
    out=ap.parse_args().out
    if out.exists():raise FileExistsError("Refuse to overwrite source audit results")
    z=compare()
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_REAL_TWO_INDEPENDENT_REPRO_STAGE_SOURCES",json.dumps({
       "Li":z["source_Li2026_reproductive_stages"]["observed_array_vs_dark_seed_per_initial_bud_stage_accounting"],
       "Boom":{k:{"stage":v["stage_accounting"],
                   "original_estimated_seed_per_flower":v["source_light_arm"]["published_original_estimated_seeds_per_flower"]}
               for k,v in z["source_Boom2020_reproductive_stages"]["source_treatment_contrasts"].items()},
       "no_grand_average":True},ensure_ascii=False))
if __name__=="__main__":main()
