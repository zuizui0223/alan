import csv
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "database/alan_experiments/arthropod_expansion_20261122/manfrin_crossendpoint"
sys.path.insert(0, str(SOURCE))
from compute_published_cpue import calculate


class TestManfrinSourceCohort(unittest.TestCase):
    def test_one_programme_not_two_studies(self):
        with (SOURCE / "programme_publication_ledger.csv").open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len({r["programme_id"] for r in rows}), 1)
        self.assertEqual(len({r["doi"] for r in rows}), 2)
        self.assertEqual(len(rows), 3)
        self.assertEqual({r["existing_physical_lighting_installation"] for r in rows}, {"same_lit_and_dark_sites"})
        self.assertEqual({r["analysis_admissibility"] for r in rows}, {"descriptive_only_unreplicated_site"})

    def test_ratio_reproduction(self):
        obj = calculate()
        self.assertAlmostEqual(obj["ratio_by_endpoint"]["aquatic_diptera_emergence"]["lit_dark_ratio_of_published_site_means"], 3.72 / 7.84)
        self.assertAlmostEqual(obj["ratio_by_endpoint"]["aquatic_diptera_flying_adults"]["lit_dark_ratio_of_published_site_means"], 416.48 / 52.17)
        self.assertAlmostEqual(obj["cross_endpoint_descriptive_ratio_of_ratios"], (416.48 / 52.17) / (3.72 / 7.84))
        self.assertLess(obj["ratio_by_endpoint"]["aquatic_diptera_emergence"]["log_ratio_descriptive_only"], 0)
        self.assertGreater(obj["ratio_by_endpoint"]["aquatic_diptera_flying_adults"]["log_ratio_descriptive_only"], 0)

    def test_physical_unit_gate(self):
        obj = calculate()
        self.assertEqual(obj["physical_exposure_sites"], {"lighted": 1, "dark": 1})
        self.assertEqual(obj["physically_independent_effects_added_to_meta"], 0)
        self.assertFalse(obj["cross_endpoint_statistical_interaction_estimated"])
        self.assertFalse(obj["within_programme_distance_to_water_treatment_interaction_estimated"])
        self.assertEqual(obj["ratio_by_endpoint"]["aquatic_diptera_emergence"]["temporal_events_not_independent_light_replicates"], 12)
        self.assertEqual(obj["ratio_by_endpoint"]["aquatic_diptera_flying_adults"]["temporal_events_not_independent_light_replicates"], 7)
        self.assertIsNone(obj["ratio_by_endpoint"]["aquatic_diptera_emergence"]["sampling_variance_for_meta"])
        self.assertIsNone(obj["ratio_by_endpoint"]["aquatic_diptera_flying_adults"]["sampling_variance_for_meta"])

    def test_frozen_output_matches(self):
        obj = calculate()
        expected = json.loads((SOURCE / "published_cpue_ratios.json").read_text(encoding="utf-8"))
        self.assertEqual(obj, expected)

if __name__ == "__main__":
    unittest.main()
