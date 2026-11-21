#!/usr/bin/env python3
"""Source-grain, non-pooling synthesis of ALAN predator and mutualist-antagonist studies.

Never treat multiple endpoints, leaves, nights, individual insects, or model
caterpillars from the same assigned light plots as independent ALAN studies.
No source whose original cluster-level variance is unavailable receives an SE.
"""
import argparse
import math
import json
from collections import Counter, defaultdict
from pathlib import Path


def synthesize(registry):
    if registry.get("schema") != "alan_20261118_new_predator_roles_programme_ledger_v1":
        raise ValueError("Unexpected frozen primary-source programme ledger schema")
    rows = registry["rows"]
    keys = [r["endpoint_id"] for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicated endpoint_id")
    if len({r["programme_id"] for r in rows}) >= len(rows):
        raise ValueError("Missing same-programme endpoint grouping")
    if any(r["source_url"].startswith("https://doi.org/") is False
           and "github.com/zuizui0223/alan" not in r["source_url"] for r in rows):
        raise ValueError("Missing verified primary source URL")
    by_ep = defaultdict(list)
    for r in rows:
        by_ep[r["endpoint"]].append(r)
        if r.get("plant_fitness_measured") not in (False, True):
            raise ValueError("Explicitly distinguish plant-fitness endpoint")
        if not r.get("source_treatment_units") or not r.get("source_treatment_units"):
            raise ValueError("No treatment unit")
        if not r.get("raw_source_status"):
            raise ValueError("No source fidelity status")
        if r["effect"].get("kind") == "published_approximate_fold_decline" and isinstance(r["effect"]["value"], float):
            raise ValueError("Cannot turn an approximate source phrase into a precise biological effect")
    predator = by_ep["predator_abundance"]
    programmes = [r["programme_id"] for r in predator]
    if len(programmes) != len(set(programmes)):
        raise ValueError("Predator direction across studies would double-count one original programme")
    dirs = Counter(r["direction"] for r in predator)
    if set(dirs) != {"increase", "decrease", "not_resolved"}:
        raise ValueError("Keep positive, negative and unresolved original predator programmes")
    # Brown (2023) published Table 3; these two within-arm summaries are NOT
    # a time-controlled DID, a between-site treatment ratio or density.
    b = next(r for r in predator if r["programme_id"] == "Brown2023")["effect"]
    if b["kind"] != "three_period_experimental_trap_means" or min(b["pre"], b["during"], b["post"]) <= 0:
        raise ValueError("Published Brown periods invalid")
    brown = {
        "experimental_trap_mean_pre": b["pre"],
        "experimental_trap_mean_during": b["during"],
        "experimental_trap_mean_post": b["post"],
        "during_to_pre_descriptive_ratio": b["during"] / b["pre"],
        "post_to_pre_descriptive_ratio": b["post"] / b["pre"],
        "during_minus_pre_insects_per_trap": b["during"] - b["pre"],
        "original_published_pre_minus_during_tukey_CI": b["reported_tukey_pre_vs_during_ci_for_pre_minus_during"],
        "not_a_control_adjusted_BACI_effect": True,
        "repeated_15_nights_are_not_independent_light_replicates": True,
    }
    # Brown's predator and herbivore guilds were classified from the SAME
    # ten experimental pitfall traps across the SAME three light periods.
    # Compute within-arm stage composition only; no predator->herbivore
    # mortality or control-adjusted causal mediation is observed.
    hrows = [r for r in by_ep["herbivore_abundance"] if r["programme_id"] == "Brown2023"]
    if len(hrows) != 1:
        raise ValueError("Brown herbivore endpoint must remain one non-independent stage")
    h = hrows[0]["effect"]
    if h["kind"] != "three_period_experimental_trap_means" or min(h["pre"], h["during"], h["post"]) <= 0:
        raise ValueError("Invalid original Brown Table 3 herbivore means")
    trophic_ratio = {stage: b[stage]/h[stage] for stage in ("pre", "during", "post")}
    brown_trophic = {
        "source": "Brown 2023 original Table 3, experimental traps only",
        "predator_average_capture": {stage: b[stage] for stage in ("pre", "during", "post")},
        "herbivore_average_capture": {stage: h[stage] for stage in ("pre", "during", "post")},
        "predator_to_herbivore_capture_ratio": trophic_ratio,
        "light_on_to_pre_predator_herbivore_ratio_of_ratios": trophic_ratio["during"]/trophic_ratio["pre"],
        "post_to_pre_predator_herbivore_ratio_of_ratios": trophic_ratio["post"]/trophic_ratio["pre"],
        "herbivore_during_to_pre_within_lit": h["during"]/h["pre"],
        "herbivore_pre_minus_during_tukey_CI_from_original": h["reported_tukey_pre_vs_during_ci_for_pre_minus_during"],
        "same_original_trap_level_process_not_two_experiments": True,
        "pitfall_counts_are_activity_density_not_true_abundance": True,
        "not_a_direct_predation_or_plant_damage_effect": True,
        "not_a_control_adjusted_BACI_estimate": True,
    }
    # Distinct endpoint: surrogate prey risk, NOT predator abundance.
    # Natural streetlights are already installed; dusk coincides with lamp
    # switching on. The original model is conditional on random lamp/batch.
    # No extra lamp clusters or unreported covariance matrix are invented.
    eck = [r for r in by_ep["surrogate_bait_predation_probability"]
           if r["programme_id"] == "Eckhartt2022"]
    if len(eck) != 1:
        raise ValueError("Eckhartt original published lamppost design missing")
    e = eck[0]["effect"]
    if e["kind"] != "published_random_lamp_batch_GLMM" or e["published_baits_total"] != 830:
        raise ValueError("Eckhartt original Table 1 model/source sample drift")
    if e["published_before_baits_predated"] + e["published_after_baits_predated"] != 516:
        raise ValueError("Eckhartt original 516 predations changed")
    if e["interaction_SE"] <= 0 or abs(e["pre_dusk_distance_interaction"]/e["interaction_SE"]-e["interaction_z"]) > .02:
        raise ValueError("Eckhartt published z and SE do not reconcile")
    ek = {
        "original_study_programme":"Eckhartt2022",
        "original_n_lamppost_clusters":7,
        "original_n_bait_events":e["published_baits_total"],
        "original_n_days":13,
        "published_model":e,
        "post_dusk_distance_logit_slope_per_m": e["distance_post_dusk_logit_per_m"],
        "pre_dusk_distance_logit_slope_per_m": e["distance_post_dusk_logit_per_m"] + e["pre_dusk_distance_interaction"],
        "post_dusk_distance_odds_ratio_per_m": math.exp(e["distance_post_dusk_logit_per_m"]),
        "pre_dusk_distance_odds_ratio_per_m": math.exp(e["distance_post_dusk_logit_per_m"] + e["pre_dusk_distance_interaction"]),
        "interaction_difference_pre_minus_post_odds_ratio_per_m": math.exp(e["pre_dusk_distance_interaction"]),
        "published_se_only_conditional_model":True,
        "no_independent_lamppost_refit":True,
        "not_randomly_assigned_ALAN":True,
        "pre_dusk_daytime_vs_post_dusk_confounded":True,
        "neither_live_caterpillar_nor_plant_fitness_measured":True,
        "cannot_meta_pool_with_predator_abundance":True,
    }
    boom = next(r for r in by_ep["nursery_seed_predation"] if r["programme_id"] == "Boom2020")["effect"]
    specs = {}
    for treatment in ("dark", "green", "white", "red"):
        a = boom[treatment]
        for v in a.values():
            if not 0 <= v <= 1:
                raise ValueError("Source proportions outside [0,1]")
        specs[treatment] = a["fruit"] * (1 - a["infested"]) * a["fertilized_noninfested"]
    if boom["seed_output_omnibus_p"] < 0.05:
        raise ValueError("Source overall seed test changed: reconsider scientific claim")
    return {
        "schema": "alan_20261118_predator_role_nonpooled_stage_synthesis_v1",
        "evidence": {
            "independent_programmes_unique": len(set(r["programme_id"] for r in rows)),
            "distinct_endpoint_entries": len(rows),
            "predator_abundance_programmes": len(predator),
            "predator_abundance_published_directions": dict(dirs),
            "opposite_predator_directions_exist": bool(dirs["increase"] and dirs["decrease"]),
            "forest_predator_direction_crosses_increase_null_decrease": set(r["direction"] for r in predator if "forest" in r["biome"]) == {"increase", "not_resolved", "decrease"},
            "no_forest_vs_open_habitat_causal_moderator_identified": True,
            "no_plant_pollination_mediation_fitted": True,
            "no_common_endpoint_random_effects_meta_fitted": True,
            "different_ecological_stages_are_not_interchangeable": sorted(by_ep),
        },
        "eckhartt2022_original_published_model": ek,
        "brown2023_original_published_period_arithmetic": brown,
        "brown2023_predator_herbivore_source_trophic_capture_reweighting": brown_trophic,
        "boom2020_spectral_nursery_algebra": {
            "product_of_published_means_fruit_times_noninfested_times_fertilized": specs,
            "green_to_dark_product_ratio": specs["green"]/specs["dark"],
            "green_to_dark_fruit_factor": boom["green"]["fruit"]/boom["dark"]["fruit"],
            "green_to_dark_noninfested_factor": (1-boom["green"]["infested"])/(1-boom["dark"]["infested"]),
            "green_to_dark_conditional_fertilization_factor": boom["green"]["fertilized_noninfested"]/boom["dark"]["fertilized_noninfested"],
            "overall_seed_output_reported_omnibus_p": boom["seed_output_omnibus_p"],
            "overall_net_seed_decline_statistically_resolved": False,
            "algebra_not_causal_mediation": True,
            "nested_ovules_not_independent_lighting_installations": True,
        },
        "programme_source_status": [{
            "programme_id": k, "endpoints": sorted(r["endpoint_id"] for r in rows if r["programme_id"] == k),
            "has_original_reanalysis_in_main": k == "Parkinson2020",
        } for k in sorted({r["programme_id"] for r in rows})],
        "scope": {
            "systematic_review": False,
            "effect_counts_not_vote_count_meta": True,
            "cannot_identify_predator_beneficial_pollination": True,
            "no_novel_experimental_results_created": True,
            "new_edi_raw_reanalysis_requires_source_verified_pair_IDs": True,
            "source_unavailable_is_not_negative_biology": True,
        },
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--registry", type=Path, default=Path("results/all_species_forest/biotic_interactions/predator_role_original_published_programme_ledger_20261118.json"))
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    ans = synthesize(json.loads(args.registry.read_text()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(ans, ensure_ascii=False, indent=2) + "\n")
    print("PREDATOR_PROGRAMMES",ans["evidence"]["independent_programmes_unique"])
    print("PREDATOR_DIRECTIONS",ans["evidence"]["predator_abundance_published_directions"])
    print("BROWN_DESCRIPTIVE_DURING_PRE",round(ans["brown2023_original_published_period_arithmetic"]["during_to_pre_descriptive_ratio"],6))
    print("BROWN_TROPHIC_PREDATOR_HERBIVORE_RATIO_OF_RATIOS",round(ans["brown2023_predator_herbivore_source_trophic_capture_reweighting"]["light_on_to_pre_predator_herbivore_ratio_of_ratios"],6))
    print("ECKHARTT_ORIGINAL_POST_DUSK_PREDATION_OR_PER_M",round(ans["eckhartt2022_original_published_model"]["post_dusk_distance_odds_ratio_per_m"],6))
    print("BOOM_SPECTRAL_PRODUCT_GREEN_DARK",round(ans["boom2020_spectral_nursery_algebra"]["green_to_dark_product_ratio"],6))


if __name__ == "__main__":
    main()
