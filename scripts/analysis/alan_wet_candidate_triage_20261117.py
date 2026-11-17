#!/usr/bin/env python3
"""Source-pinned candidate-mechanism *triage*, not causal mediation or a pooled meta-effect.

Run at repository root:
    python scripts/analysis/alan_wet_candidate_triage_20261117.py

All inputs are existing alan/main aggregate tables, never reselected using signs.
Python stdlib only. A report with existing results is not a new independent study.
"""
import csv
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results" / "all_species_forest"


def rows(path):
    with (BASE / path).open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file))


def find_one(data, **selectors):
    out = [r for r in data if all(r.get(k) == v for k, v in selectors.items())]
    if len(out) != 1:
        raise AssertionError(f"Expected one row for {selectors}, found {len(out)}")
    return out[0]


def pearson(x, y):
    xbar, ybar = sum(x) / len(x), sum(y) / len(y)
    cross = sum((a-xbar)*(b-ybar) for a,b in zip(x,y))
    ssx = sum((a-xbar)**2 for a in x)
    ssy = sum((b-ybar)**2 for b in y)
    return cross / math.sqrt(ssx * ssy)


def ranks(v):
    # No ties among the four original species with independently known EIV.
    if len(set(v)) != len(v):
        raise AssertionError("Rank tie requires mid-rank treatment")
    return [sorted(v).index(z)+1 for z in v]


def inspect():
    macro = rows("wet_environment/environment.csv")
    regional = rows("wet_environment/regions.csv")
    species = rows("wet_environment/giavi_diurnal.csv")
    traits = rows("wet_drivers/traits.csv")
    interactions = rows("wet_drivers/site_interactions.csv")
    prediction = rows("wet_drivers/regional_prediction.csv")
    annual = rows("wet_annual/contrasts.csv")
    frozen = find_one(macro, model="frozen", term="candidate")
    adjusted = find_one(macro, model="environment", term="candidate")
    # Halt if underlying source cohort/definition silently changes.
    assert (int(adjusted["candidates"]), int(adjusted["references"]),
            int(adjusted["regions"])) == (158, 1419, 12)
    assert abs(float(frozen["estimate"]) - 0.5922104455502863) < 1e-8
    assert abs(float(adjusted["estimate"]) - 0.5338399813053877) < 1e-8

    positive = [r for r in regional if float(r["difference"]) > 0]
    negative = [r for r in regional if float(r["difference"]) < 0]
    assert len(regional) == 12 and len(positive) == 7 and len(negative) == 5
    assert sum(int(r["candidates"]) for r in regional) == 158
    ranked = sorted(regional, key=lambda r: int(r["candidates"]), reverse=True)
    top2 = sum(int(r["candidates"]) for r in ranked[:2])
    big = [r for r in regional if int(r["candidates"]) >= 10]

    measured = [r for r in species if r["moisture"].strip() and r["estimate"].strip()]
    assert len(species) == 21 and len(measured) == 4
    x = [float(r["moisture"]) for r in measured]
    y = [float(r["estimate"]) for r in measured]
    xy = list(zip(x,y))
    wet = [b for a,b in xy if a > 7]
    other = [b for a,b in xy if a <= 7]
    assert len(wet) == 2 and len(other) == 2
    sx, sy = ranks(x), ranks(y)
    rho = pearson(sx,sy)
    exact_extreme = sum(
        abs(pearson(sx, list(p))) >= abs(rho)-1e-12
        for p in itertools.permutations(sy)
    )
    slope = sum((a-sum(x)/len(x))*(b-sum(y)/len(y)) for a,b in xy) / (
        sum((a-sum(x)/len(x))**2 for a in x)
    )
    strict = [r for r in measured if int(r["shared_meadow_pairs"]) >= 2]
    strict_wet = [r for r in strict if float(r["moisture"]) > 7]
    strict_other = [r for r in strict if float(r["moisture"]) <= 7]
    assert len(strict) == 2 and len(strict_wet) == 2 and not strict_other

    t_q = min(float(r["q"]) for r in traits)
    all_site = [r for r in interactions if r["scope"] == "all"]
    b_site = [r for r in interactions if r["scope"] == "brightening"]
    assert len(traits) == 34 and len(all_site) == 12 and len(b_site) == 12
    baseline = find_one(prediction, omitted="all", model="baseline")
    added = find_one(prediction, omitted="all", model="plus_traits")
    annual_effect = find_one(annual, scope="consecutive_surveys",
                             state="abs_water_freq", outcome="difference")
    near = find_one(interactions, scope="all", modifier="near")

    return {
        "schema": "alan_wet_mechanism_candidate_triage_v1",
        "scope": "secondary post-hoc mechanism prioritization; no new primary ALAN effect",
        "primary_existing_moisture_difference": {
            "n_candidate": 158, "n_reference": 1419, "regions": 12,
            "unadjusted": float(frozen["estimate"]),
            "adjusted": float(adjusted["estimate"]),
            "adjusted_interval": [float(adjusted["low"]), float(adjusted["high"])]
        },
        "region_transport_diagnostic": {
            "positive": len(positive), "negative": len(negative),
            "top_two_regions_candidate_share": top2 / 158,
            "regions_at_least_10_candidates": len(big),
            "positive_among_10plus_regions": sum(float(r["difference"]) > 0 for r in big),
            "note": "same selected species; not independent geographic replication"
        },
        "previous_macro_falsifiers": {
            "extra_trait_count": len(traits), "best_adjusted_q_extra_traits": t_q,
            "regional_AUC_baseline": float(baseline["within_region_auc"]),
            "regional_AUC_plus_34_traits": float(added["within_region_auc"]),
            "site_modifiers_count": len(all_site),
            "best_site_modifier_q": min(float(r["q"]) for r in all_site),
            "brightening_only_site_modifier_count": len(b_site),
            "best_brightening_modifier_q": min(float(r["q"]) for r in b_site),
            "near_water_interaction": float(near["estimate"]),
            "near_water_interaction_q": float(near["q"]),
            "wetland_annual_study_units": int(annual_effect["sites"]),
            "wetland_consecutive_year_effect": float(annual_effect["estimate"]),
            "wetland_consecutive_year_interval": [float(annual_effect["low"]), float(annual_effect["high"])]
        },
        "swiss_same_experiment_moisture_pilot": {
            "experiment_count": 1,
            "all_plant_taxa": len(species),
            "wet_indicator_known_species": len(measured),
            "known_species": [
                {"species": r["label"].replace("_", " "),
                 "moisture": float(r["moisture"]),
                 "light_day_visit_log1p_effect": float(r["estimate"]),
                 "shared_meadow_pairs": int(r["shared_meadow_pairs"])}
                for r in measured
            ],
            "continuous_ols_per_moisture_point": slope,
            "pearson_r_descriptive": pearson(x,y),
            "spearman_r_descriptive": rho,
            "spearman_exact_permutation_two_sided": exact_extreme/math.factorial(len(sx)),
            "wet_gt7_minus_drier_species_mean_effect": sum(wet)/2-sum(other)/2,
            "species_with_at_least_two_meadow_pairs": len(strict),
            "drier_species_with_at_least_two_pairs": len(strict_other),
            "eligible_as_cross_species_effect_meta": False,
            "note": "4 species, 1 illumination experiment, species estimates share meadow installations; 0 drier taxa pass 2-pair gate. No causal generalization"
        },
        "meta_readiness": {
            "homogeneous_multi_study_lit_dark_plant_response_with_wet_dry_moderator": False,
            "effect_size": None,
            "reason": "Missing matched wet/dry taxon outcomes and sampling variances from independent illuminated experiments"
        }
    }


def main():
    result = inspect()
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))
    return result


if __name__ == "__main__":
    main()
