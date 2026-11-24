#!/usr/bin/env python3
"""Audit descriptive ratios of *published* Manfrin 2025 site-summed CPUE.

These are published group mean/SD, NOT original experimental units. No meta SE
or light-induced causal population effect may be derived from temporal samples.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "published_cpue.csv"
DEST = ROOT / "published_cpue_ratios.json"
PROGRAMME_ID = "Westhavelland_2012_2014"


def calculate():
    with SRC.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 2 or {r["endpoint"] for r in rows} != {
        "aquatic_diptera_emergence", "aquatic_diptera_flying_adults"
    }:
        raise ValueError("expected exactly two different source endpoints")
    result = {}
    for row in rows:
        if row["programme_id"] != PROGRAMME_ID or row["calendar_year"] != "2014":
            raise ValueError("wrong programme or observational year")
        if (int(row["lighted_site_count"]), int(row["dark_site_count"])) != (1, 1):
            raise ValueError("this audit does not cover replicated sites")
        if row["meta_eligible"].lower() != "false":
            raise ValueError("single treated and control site cannot establish independent-lamp meta variance")
        a, b = float(row["lit_mean_cpue"]), float(row["dark_mean_cpue"])
        if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
            raise ValueError("non-positive or invalid published means")
        result[row["endpoint"]] = {
            "lit_dark_ratio_of_published_site_means": a / b,
            "log_ratio_descriptive_only": math.log(a / b),
            "lit_mean": a,
            "lit_sd_temporal_observation": float(row["lit_sd_cpue"]),
            "dark_mean": b,
            "dark_sd_temporal_observation": float(row["dark_sd_cpue"]),
            "temporal_events_not_independent_light_replicates": int(row["temporal_sampling_events"]),
            "sampling_variance_for_meta": None,
        }
    a = result["aquatic_diptera_emergence"]["lit_dark_ratio_of_published_site_means"]
    b = result["aquatic_diptera_flying_adults"]["lit_dark_ratio_of_published_site_means"]
    return {
        "programme_id": PROGRAMME_ID,
        "source_article_doi": "10.1007/s00027-025-01161-7",
        "source_basis": "published_site_aggregated_means_and_temporal_sd_only",
        "physical_exposure_sites": {"lighted": 1, "dark": 1},
        "ratio_by_endpoint": result,
        "cross_endpoint_descriptive_ratio_of_ratios": b / a,
        "cross_endpoint_descriptive_log_ratio_difference": math.log(b) - math.log(a),
        "cross_endpoint_statistical_interaction_estimated": False,
        "within_programme_distance_to_water_treatment_interaction_estimated": False,
        "physically_independent_effects_added_to_meta": 0,
        "warnings": [
            "two different CPUE denominators and different temporal sampling schedules",
            "temporal SDs and 12/7 sampling events cannot estimate the light-installation sampling variance",
            "2017 whole-aquatic-insect emergence and 2014 Diptera-only emergence are distinct year/taxon endpoints",
            "one treatment site versus one dark site makes site and light perfectly confounded",
            "capture activity can rise without growth of the source population",
            "no pollen or wet-plant fitness response was co-measured",
        ],
    }


def main():
    obj = calculate()
    DEST.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("Manfrin source panel: 1 programme; 2 descriptive endpoints; 0 meta-eligible effects")
    print("Emergence ratio:", obj["ratio_by_endpoint"]["aquatic_diptera_emergence"]["lit_dark_ratio_of_published_site_means"])
    print("Flying ratio:", obj["ratio_by_endpoint"]["aquatic_diptera_flying_adults"]["lit_dark_ratio_of_published_site_means"])

if __name__ == "__main__":
    main()
