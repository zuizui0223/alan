#!/usr/bin/env python3
"""Audit independently confirmed source chronology mismatch (Wonderlin study).

Original author's study-year tag is 2022, while published Wonderlin
(Ecological Entomology, 2024, DOI 10.1111/een.13296) states its floral visitor
videos were summer 2022 and its actual day/night reproductive bagging trial was
summer 2019. This is *source metadata*, not a new site/year-light assignment.

Strictly refit the unchanged 301-source historical-light night reproductive
association after excluding the ONE source publication with this documented
year discordance. No causal light inference and no new corrected PCNL values
are claimed; zero individual plant/site/year records are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, estimate, require, LIGHT
)

FLAGGED_STUDY = "Wonderlin_2023"
PUBLISHED_ARTICLE = "10.1111/een.13296"
RECORDED_AUTHOR_SOURCE_YEAR = 2022
PUBLISHED_ACTUAL_BAGGING_YEAR = 2019
PUBLISHED_VIDEO_OBSERVATION_YEAR = 2022

FORMULA = ("night_share ~ " + LIGHT +
           " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")


def condensed(obj):
    return {
        "beta": obj.get("pcnl_prior_light_coefficient_fraction_per_mean_log1p"),
        "CI95_study_cluster": obj.get("95CI_study_cluster"),
        "p_nominal": obj.get("p_nominal_study_cluster"),
        "pairs": obj.get("n_paired_original_experiment_locations_year_metrics"),
        "independent_studies": obj.get("n_independent_original_studies"),
        "species": obj.get("n_source_plant_species"),
        "study_species_two_way": obj.get("two_way_independent_study_plant_cluster"),
        "status": obj.get("status"),
    }


def audit(joined, receipt, out):
    p = Path(out)
    require(not p.exists(), "Existing source audit may not be overwritten")
    d, source = load_source(joined, receipt)
    require((len(d), d.study_key.nunique(), d.species.nunique()) ==
            (301, 118, 123), "Original global field reproductive cohort changed")
    w = d.loc[d.study_key.eq(FLAGGED_STUDY)]
    require(len(w) >= 1, "Original source misdated study not found")
    require(w.year_exp.eq(RECORDED_AUTHOR_SOURCE_YEAR).all(),
            "Original source year no longer 2022")
    baseline = estimate(d, "original_full_year_model", FORMULA)
    require(baseline["status"].startswith("EXPLORATORY_FIELD_STUDY"),
            "Unchanged original source slope not estimable")
    require(abs(baseline["pcnl_prior_light_coefficient_fraction_per_mean_log1p"]
                - .0030676) < 5e-6,
            "Baseline prior-PCNL source slope does not match original frozen result")
    clean = d.loc[~d.study_key.eq(FLAGGED_STUDY)].copy()
    filtered = estimate(clean, "exclude_chronology_discordant_paper", FORMULA)
    require(filtered["status"].startswith("EXPLORATORY_FIELD_STUDY"),
            "Strict time-diagnostic omission did not estimate")
    output = {
        "schema": "ALAN_original301_publication_verified_reproductive_year_mismatch_sensitivity_v1",
        "status": "REAL_SOURCE_301_CHRONOLOGY_DISCORDANT_STUDY_OMISSION_REFIT",
        "provenance": {
            "original_preexperiment_LiCao_annual_PCNL_years": [1992, 2021],
            "annual_raster_MD5s_verified": source["annual_rasters_MD5_source_verified"],
            "original_global_field_source_SHA256": source["joined_file_sha256"],
            "source_full_original_DvN_paired": 301,
            "same_experiment_reproductive_outcome_metadata_not_relabelled": True,
        },
        "case": {
            "original_published_study_reference": FLAGGED_STUDY,
            "original_author_year_field": RECORDED_AUTHOR_SOURCE_YEAR,
            "independent_published_article_DOI": PUBLISHED_ARTICLE,
            "published_bagging_experiment_year": PUBLISHED_ACTUAL_BAGGING_YEAR,
            "published_later_video_visitor_year": PUBLISHED_VIDEO_OBSERVATION_YEAR,
            "published_video_and_experiment_are_NOT_same_year": True,
            "original_year_could_include_post_experiment_light_if_used_as_start": True,
            "original_contrasts_affected_by_study_exclusion": int(len(w)),
            "number_independent_source_studies_removed": 1,
            "we_did_not_impute_corrected_pcNL_without_original_annual_pixel_matching": True,
        },
        "study_equal_adjusted_historical_light_vs_actual_reproduction": {
            "original301": condensed(baseline),
            "excluding_independently_verified_year_discordant_paper": condensed(filtered),
            "full_minus_omission_beta": float(
                baseline["pcnl_prior_light_coefficient_fraction_per_mean_log1p"]-
                filtered["pcnl_prior_light_coefficient_fraction_per_mean_log1p"])
        },
        "guardrails": {
            "original_study_source_2022_may_be_visitor_2022_not_trial_2019": True,
            "no_reading_of_narrative_as_original_pollen_deposition": True,
            "no_assumption_87_DOI_metadata_links_identify_field_visitors": True,
            "comparison_not_a_causal_effect_of_ALAN": True,
            "preexisting_PCNLPriors_restored_only_no_raster_redownload": True,
            "all_original_individual_species_site_coordinates_withheld": True,
        },
    }
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(output, indent=2, ensure_ascii=False,
                            allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_EXPERIMENT_YEAR_WONDERLIN_SENSITIVITY",
          json.dumps(output, ensure_ascii=False), flush=True)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("joined", "receipt", "out"):
        parser.add_argument("--"+name, required=True, type=Path)
    a = parser.parse_args()
    audit(a.joined, a.receipt, a.out)
