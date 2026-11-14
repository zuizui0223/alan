#!/usr/bin/env python3
"""Cross-source, evidence-gated ALAN vulnerability research audit.

Mechanistic objective: ALAN × nocturnal *effective* pollination dependence ×
selfing/daytime functional rescue × predator's negative or positive function ×
flower/pollinator phenological overlap × island/mainland, NOT individual
species-level GloBI labels.

Uses independently archived public aggregate summaries, signed-source SHA pins
and a public article metadata registry. No undocumented local guild is coded 0.
This audit produces a *feasibility and inferential guardrail* report, not a
new unsupported many-way regression or a synthetic meta-analysis.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

KNOWN_SOURCE_SHA = "670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3"
KNOWN_PHOTO_SHA = "f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c"
KNOWN_INS_SHA = "64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852"


def require(ok, reason):
    if not bool(ok):
        raise ValueError("ALAN_EVIDENCE_GATE: " + reason)


def read_json(root, rel):
    p = Path(root)/rel
    require(p.is_file(), f"Missing original documented aggregate {rel}")
    return json.loads(p.read_text(encoding="utf-8"))


def contains(root, rel, literals):
    path=Path(root)/rel
    require(path.exists(), f"Missing independently verified original source note {rel}")
    s=path.read_text(encoding="utf-8")
    for v in literals:
        require(v in s, "Expected original-source note changed: "+v)


def finite_ci(obj, key):
    x=obj[key]
    require(isinstance(x,list) and len(x)==2 and
            all(isinstance(a,(int,float)) and math.isfinite(a) for a in x) and
            x[0]<=x[1], "Invalid original uncertainty bound: "+key)
    return list(x)


def main(root, out):
    out=Path(out)
    require(not out.exists(), "Never overwrite the evidence receipt")
    reg=read_json(root,"config/alan_night_syndrome_multichannel_macro_sources_20261114.json")
    require(reg["schema"]=="ALAN_source_locked_macro_syndrome_predation_phenology_island_registry_v1","Article registry schema drift")
    require(len(reg["hypotheses"])==4,"Four scientific domains changed")
    require(all(reg["rules"].values()),"At least one source false absence shortcut allowed")
    result={}
    gb=read_json(root,"results/alan_v011_moth_butterfly_flower_documentation_multi_channel_ALAN_20261112.json")
    pl=read_json(root,"results/alan_globi_v011_pl64_readme_moderator_robustness_20261111.json")
    pred=read_json(root,"results/alan_GloBI_v011_predator_dual_role_evidence_tiers_20261111.json")
    rept=read_json(root,"results/alan_true13yr_moth_gecko_sameisland_light_source_verified_20261112.json")
    long=read_json(root,"results/alan_GloBI_long_photo_2012_2025_source_verified_20261111.json")
    morph=read_json(root,"results/alan_realGBIF_gecko_moth_photo_light_spatial_temporal_parallel_20261112.json")
    require(gb["source"]["full_GloBI_v011_SHA256"]==KNOWN_SOURCE_SHA and
            pred["source"]["original_interactions_sha256"]==KNOWN_SOURCE_SHA and
            gb["source"]["original_photo_SHA256"]==KNOWN_PHOTO_SHA and
            long["source"]["historical_original_photo_SHA256"]==KNOWN_PHOTO_SHA,
            "Global original source SHA identity mismatch")
    require(rept["source"]["insect_source_joined_csv_sha256"]==KNOWN_INS_SHA,
            "Animal original GBIF insect source SHA changed")
    require(gb["source"]["photo_site_years"]==10698 and
            gb["source"]["photo_island_site_years"]==933 and
            gb["source"]["photo_mainland_site_years"]==9765,
            "Original mainland-island source opportunity counts changed")
    require(pred["source"]["full_interaction_rows"]==24453954,
            "GloBI complete original source count drift")
    earlier=pl["exact_old_model_reproduction"]
    require(abs(earlier["original_README_PL64_beta"]-(-.01938004575601668))<1e-12 and
            abs(earlier["original_README_PL64_p"]-.043993409104236464)<1e-12 and
            abs(gb["source"]["original_legacy_ALAN_slope_reconstructed"]-
                earlier["moth_positive_legacy_v0_broad_beta"])<1e-9,
            "Frozen original positive source estimates do not match")
    test30=next(x for x in pl["actual_model_estimates"] if
                x["cohort"]=="GloBI v0.11 moth positive")
    require(test30["n_species"]==30 and abs(test30["p_nominal"]-.4799746115289204)<1e-10,
            "Original v0.11 GloBI source reanalysis changed")
    # Detect a real metadata-only old field that accidentally contains PL30's
    # p value in a displayed README PL64 slot. Do NOT silently rewrite source.
    earlier_labels=gb["frozen_positive_result_unmodified"]
    wrong_repeated_p=(
        abs(earlier_labels["p_nominal"]-test30["p_nominal"])<1e-10 and
        abs(earlier_labels["p_nominal"]-earlier["original_README_PL64_p"])>.01
    )
    require(wrong_repeated_p,"Historical label metadata bug no longer exists; refresh cross-source provenance handling")

    contains(root,"docs/ALAN_ORIGINAL_DAYNIGHT_REPRODUCTIVE_562ISLAND_PRIORLIGHT_20261113.md",
             ["**562 exact GSHHG islands", "**52** of 134", "**2,373 cell", "39 independent islands"])
    contains(root,"docs/ALAN_GLOBAL562_DIRECT_FLORAL_MOTHLIKE_PCNL_EMPIRICAL_20261113.md",
             ["14,330 distinct species", "**1,691 species**", "**257 species**", "**13 species**",
              "genus fixed effects"])
    contains(root,"docs/ALAN_301_ORIGINAL_VISITOR_OBS_OA_SAMPLING_SELECTION_20261114.md",
             ["**92**", "247 of the 301", "**36 candidate OA links**"])
    contains(root,"docs/ALAN_TRUE_WENZELL2019_POPULATIONYEAR_GUILD_REPRODUCTION_LINK_20261114.md",
             ["SIX of these seven", "5/6 mapped", "within **1 km**"])
    contains(root,"docs/ALAN_ORIGINAL_KAYS2012_GAMBOA_KINKAJOU_FRUIT_CLOSED_SEMANTICS_20261114.md",
             ["**bees can pass**","2010", "same"])

    study_sources={
       "island_macro_observed_photo_opportunity":{"independent_islands":562,
         "photo_cell_years":77647,"accepted_recorded_species":14330,
         "measured_real_day_night_source_species_in_photo_frame":52,
         "species_with_direct_color_and_shape":1691,
         "white_tubular_like_source_species":257,
         "white_tube_deep_source_species":13,
         "original_species_site_year_observations":2012925,
         "phenology_annual_field_plant_peak_or_local_visitor_activity_verified":False,
         "predator_abundance_or_moth_capture_at_same_flower_measured":False},
       "island_mainland_photo_alanspatial":{"island_siteyears":933,
         "mainland_siteyears":9765,
         "geographical_5degree_blocks":78,
         "full_GloBI_positive_source_records":24453954,
         "moth_only_documented_not_daytime_absence":True,
         "slope_both_vs_moth_only":gb["analysis"]["direct_planned_contrasts"]["both_documented_minus_moth_only_recorded"],
         "original_legacy_moth_light_beta":earlier["moth_positive_legacy_v0_broad_beta"],
         "label_updated_moth_light_beta":0.750257336},
       "reproductive_field_night_only_share_vs_historical_light":{"pairs":301,
         "independent_papers":118,
         "experiment_study_year_local_light_randomized":False,
         "parent_pcnl_beta":.003067606366053171,
         "parent_pcnl_ci95":[-.02238802830826334,.02852324104036968],
         "parent_pcnl_p_nominal":.8132842743524381,
         "original_papers_with_some_visitor_observations":92,
         "papers_with_local_year_plant_guild_and_reproduction_jointly_confirmed_at_strict_tier":2,
         "no_ever_equate_source_301_and_562_at_site_year":True},
       "animal_photo_predator_and_moth":{"matched_same_island_cell_years":15255,
         "primary_model_islands":rept["primary"]["independent_islands"],
         "primary_cell_years":rept["primary"]["photo_cell_years"],
         "gecko_fraction_light_moth_photo_interaction_beta":rept["primary"]["gecko_fraction_by_priorloglight_interaction"]["beta"],
         "interaction_ci95":rept["primary"]["gecko_fraction_by_priorloglight_interaction"]["CI95"],
         "photographed_reptile_not_predator_attack":True},
       "predator_positive_named_roles":{
         "GloBI_strict_bats_eat_moth_and_db_pollinates_species":pred["group_counts"]["bats"]["source_restricted"],
         "GloBI_strict_geckos_eligible_species":pred["group_counts"]["geckos"]["source_restricted"],
         "GloBI_strict_centipedes_eligible_species":pred["group_counts"]["centipedes"]["source_restricted"],
         "matched_local_actual_both_roles_at_one_flower":0,
         "source_positive_not_inferred_as_effective_pollen":True},
       "phenology_observation_coverage":{"photo_observation_years":14,
         "original_photo_siteyears":long["counts"]["source_photo_site_years_with_recognized_landmass"],
         "original_photo_year_region_matched_island_mainland":long["counts"]["year_region_island_mainland_matched_source_site_years"],
         "early_2014_island_sites_with_matched_mainland_5degree_source":0,
         "source_annual_light_attached_to_pre_2019_photos":False,
         "direct_flowering_pollinator_phase_overlap_estimated":False},
    }
    require(study_sources["animal_photo_predator_and_moth"]["primary_model_islands"]==50 and
            study_sources["animal_photo_predator_and_moth"]["primary_cell_years"]==3424 and
            study_sources["phenology_observation_coverage"]["original_photo_siteyears"]==332576,
            "Historic frozen source count changed")
    finite_ci(study_sources["island_mainland_photo_alanspatial"]["slope_both_vs_moth_only"],"CI95")
    finite_ci(study_sources["animal_photo_predator_and_moth"],"interaction_ci95")
    # Distinguish denominator-limited direct function from high image opportunity.
    photo=study_sources["island_macro_observed_photo_opportunity"]
    direct_fraction=photo["measured_real_day_night_source_species_in_photo_frame"]/photo["accepted_recorded_species"]
    morph_fraction=photo["species_with_direct_color_and_shape"]/photo["accepted_recorded_species"]
    require(0<direct_fraction<.01 and .05<morph_fraction<.25,
            "Wrong primary source cohort denominators")
    # The 2026 primary experiment has real light and observed visitor change,
    # yet does not have actual night/day seed-set for all 29/40 syndrome plants.
    papers=reg["published_sources"]
    qtp=[x for x in papers if x["doi"]=="10.1002/ece3.74368"][0]
    geo=[x for x in papers if x["doi"]=="10.1093/jpe/rtag160"][0]
    swiss=[x for x in papers if x["doi"]=="10.1016/j.biocon.2025.111650"][0]
    require(qtp["independent_spatial_blocks"]==5 and
            qtp["night_lepidoptera_pct_change"]==-57.5 and
            geo["plant_species"]==2 and
            swiss["flowering_schedule_species_n"]==16,
            "Published experimental-source registry metrics were altered")
    # NO four-factor regression is identifiable from separate summary sources.
    axis=[
      {"id":"H1","status":"PARTIAL_EXPERIMENTAL_EVIDENCE_NO_GLOBAL_COMMON_SUPPORT",
       "real_manipulated_ALAN_examples":1,
       "source_measured_plant_species_in_light_reproduction_experiment":2,
       "global_562_island_observed_species_experimentally_night_documented":52,
       "reason":"Light-treated night/seed reproduction exists only in a tiny plant subset; field 301 historic light is observational and different study sites"},
      {"id":"H2","status":"NOT_ESTIMABLE_DIRECT_PREDATION_AND_COMPENSATORY_POLLEN",
       "same_flower_contemporaneous_predator_attack_and_effective_pollen_with_light":0,
       "reason":"GBIF reptiles are positive photographs; GloBI predator interactions are not adult same-site capture and pollen deposition"},
      {"id":"H3","status":"NOT_ESTIMABLE_DIRECT_SEASONAL_OR_DIEL_MISMATCH",
       "direct_complete_local_flower_pollinator_phenophase_pairs_under_light":0,
       "reason":"Swiss flower schedule manipulation and Tibetan guild network record different populations; photo date is observer timing, not flower phenophase"},
      {"id":"H4","status":"NOT_ESTIMABLE_EFFECT_OF_PRIOR_ISLAND_POLLINATOR_SHIFT",
       "verified_global_prior_shifts_connected_to_light_and_real_reproduction":0,
       "reason":"Syndrome proxy and documented global visits cannot be assumed an island population has already shifted from the ancestor"}
    ]
    require([x["id"] for x in axis]==[x["id"] for x in reg["hypotheses"]],
            "Predeclared hypothesis ordering mismatch")
    result={
      "schema":"ALAN_20261114_macro_mechanism_assessment_actual_source_gate_v1",
      "status":"REAL_SAVED_AGGREGATE_SOURCE_QA_ONLY_NO_NEW_INTERACTION_BETA",
      "original_source_aggregate_replay":study_sources,
      "evidence_support_conversion":{
          "true_direct_reproductive_night_characterized_species_fraction_of_14330":direct_fraction,
          "direct_shape_color_recorded_species_fraction_of_14330":morph_fraction,
          "source_photographic_species_does_not_equal_local_moth_dep":True,
          "prior_52_reproductive_species_and_257_morph_candidates_NOT_proven_common_species_intersection":True},
      "conditional_hypotheses_identification":axis,
      "experimental_independent_source_examples":[
        {"doi":qtp["doi"],"n_blocks":qtp["independent_spatial_blocks"],
         "night_moth_visitation_percent_change":qtp["night_lepidoptera_pct_change"],
         "night_diptera_visitation_percent_change":qtp["night_diptera_pct_change"],
         "not_measured_same_plant_mature_seed_for_all":True},
        {"doi":geo["doi"],"source_plant_species":geo["plant_species"],
         "SI_night_dependent_seed_change_percent":geo["SI_nocturnal_night_ALAN_seed_set_pct_change"],
         "same_station_correlated_with_network_trial":True},
        {"doi":swiss["doi"],"lit_sites":swiss["illuminated_sites"],
         "dark_sites":swiss["paired_controls"],
         "species_with_flowering_schedule_changes_percent":swiss["flowering_schedule_species_in_percent"],
         "actual_pollinator_mismatch_NOT_measured":True},
      ],
      "old_source_record_provenance_bug":{
         "file":"results/alan_v011_moth_butterfly_flower_documentation_multi_channel_ALAN_20261112.json",
         "stored_field":"frozen_positive_result_unmodified.p_nominal",
         "value_currently_written":earlier_labels["p_nominal"],
         "actually_readme_PL64_p_from_direct_original_replay":earlier["original_README_PL64_p"],
         "incorrectly_copied_v011_PL30_p":test30["p_nominal"],
         "prior_original_file_not_overwritten":True},
      "safe_inference":{
         "all_prior_original_frozen_source_positive_results_preserved":True,
         "not_claimed_as_meta_analytic_global_combined_ALAN_effect":True,
         "no_absent_GBIF_GloBI_trait_as_biological_absence":True,
         "no_cross_scale_beta_subtraction":True,
         "remote_sensing_light_not_local_lux":True,
         "published_field_experiment_effects_separate_from_observational_GBIF":True,
         "no_favorable_subgroup_sign_selection":True,
         "local_phenology_and_syndrome_shift_still_missing_for_macro_4way":True,
      },
      "predeclared_next_phase":{
         "unit":"matched taxon-or-functional-guild × locality (island and mainland) × sampling year × month AND diel time window",
         "exposure":"measured lamp lux/spectrum where possible; otherwise frozen historical PCNL or VIIRS remote source with separate interpretation",
         "primary_outcome":"measured fruit/seed or conspecific stigma pollen under day/night/open/closed, NOT iNat photo proportion",
         "secondary_outcome":"effort-standardized pollinator visitation/guild transfer and plant phenology synchronization",
         "repeated_geo":"physical GSHHG island; nearest mainland/lat-climate analogous blocks; same species/lineages as available",
         "predator":"adult moth prey attacks, negative fear effect, predator's directly observed pollen-transfer service, each as separate variables",
         "mismatch":"daily flowering-open and pollinator foraging time + seasonal onset/peak and activity; DOY not observation upload time",
         "tests":["ALAN×night_dependence×effective_redundancy","ALAN×predator_actual_negative_plus_positive",
                  "ALAN×phenological_overlap","ALAN×island×verified_preexisting_shift"],
         "minimum_independent_papers_per_functional_class":25,
         "minimum_unique_plant_species_per_functional_class":20,
         "report_not_estimable_until_real_data_satisfy_gate":True,
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_TRUE_MACRO_SYNDROME_SOURCE_GATES",
          json.dumps({"source":study_sources,
                      "direct_trait_coverage":result["evidence_support_conversion"],
                      "estimands":axis,
                      "old_p_transcription":result["old_source_record_provenance_bug"],
                      "status":result["status"]},ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,default=Path("."))
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    main(a.repo,a.out)
