#!/usr/bin/env python3
"""Non-invention guardrails for known published visitor-vs-treatment year mismatch."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_reproduction_published_year_audit_20261114 import (
    FLAGGED_STUDY, PUBLISHED_ARTICLE, RECORDED_AUTHOR_SOURCE_YEAR,
    PUBLISHED_ACTUAL_BAGGING_YEAR, PUBLISHED_VIDEO_OBSERVATION_YEAR, condensed
)


class OriginalFieldReproductionChronology(unittest.TestCase):
    def test_source_discordance_not_permitted_as_local_same_year_visitor(self):
        self.assertEqual(FLAGGED_STUDY, "Wonderlin_2023")
        self.assertEqual(PUBLISHED_ARTICLE, "10.1111/een.13296")
        self.assertEqual(RECORDED_AUTHOR_SOURCE_YEAR, 2022)
        self.assertEqual(PUBLISHED_ACTUAL_BAGGING_YEAR, 2019)
        self.assertEqual(PUBLISHED_VIDEO_OBSERVATION_YEAR, 2022)
        self.assertNotEqual(PUBLISHED_ACTUAL_BAGGING_YEAR,
                            PUBLISHED_VIDEO_OBSERVATION_YEAR)

    def test_result_subset_does_not_export_source_rows(self):
        x=condensed({
            "pcnl_prior_light_coefficient_fraction_per_mean_log1p": .01,
            "95CI_study_cluster": [-.1, .2],
            "p_nominal_study_cluster": .6,
            "n_paired_original_experiment_locations_year_metrics": 300,
            "n_independent_original_studies": 117,
            "n_source_plant_species": 122,
            "two_way_independent_study_plant_cluster": {"p_two_sided_t":.67},
            "status": "EXPLORATORY_FIELD_STUDY_NIGHT_VS_DAY_REPRODUCTIVE_OUTPUT_PRIOR_PCNL",
            "lat": 11.1,
            "lon": 22.2,
            "species": "Hidden_Plant"
        })
        self.assertEqual(x["pairs"],300)
        self.assertEqual(x["species"],122)
        self.assertNotIn("lat", x)
        self.assertNotIn("lon", x)


if __name__=="__main__":
    unittest.main()
