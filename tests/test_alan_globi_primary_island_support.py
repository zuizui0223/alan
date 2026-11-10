"""Tests for GloBI-first source-support logic (synthetic evidence, no true absences)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts" / "analysis"))
import alan_globi_primary_island_support as m


class TestGloBIPrimary(unittest.TestCase):
    def setUp(self):
        self.globi = pd.DataFrame([
            ("Flora alba", "moth_family_visit", "visitsFlowersOf", "research", "Sphingidae"),
            ("Flora alba", "butterfly_family_visit", "visitsFlowersOf", "research", "Nymphalidae"),
            ("Flora rubra", "moth_family_visit", "visitsFlowersOf", "iNaturalist", "Noctuidae"),
            ("Flora rosea", "butterfly_family_visit", "visitsFlowersOf", "research", "Lycaenidae"),
        ], columns=["plant_name", "guild", "relation", "study_title", "source_family"])
        self.flora = pd.DataFrame([
            ("A", "Flora alba", "introduced", "False"),
            ("A", "Flora rubra", "native", "False"),
            ("A", "Flora rosea", "unresolved", "False"),
            ("B", "Flora alba", "native", "False"),
        ], columns=["island_id", "accepted_species", "origin_status", "status_conflict"])
        self.photos = pd.DataFrame([
            ("A", "2020", "5", "3", "1", "quarter_months_NOT_annual"),
            ("B", "2021", "8", "1", "0", "NO_same_year_light"),
        ], columns=["island_id", "year", "photos_all_original_no_roi_gate",
                    "photos_with_positive_moth_visit_record",
                    "same_year_viirs_matched_photos", "light_source"])

    def test_moth_and_butterfly_can_overlap(self):
        a = m.scenario(self.globi, self.flora, self.photos)
        self.assertEqual(a["distinct_moth_positive_plant_names"], 2)
        self.assertEqual(a["plants_with_both_recorded_guilds"], 1)
        self.assertEqual(a["exact_moth_plant_name_matches_in_island_flora"], 2)
        self.assertEqual(a["moth_plant_island_name_rows"], 3)
        self.assertEqual(a["introduced_origin_moth_species_NOT_naturalization"], 1)

    def test_source_exclusion_changes_coverage_not_ecology(self):
        g = self.globi[~self.globi.study_title.str.contains("inaturalist", case=False)]
        a = m.scenario(g, self.flora, self.photos)
        self.assertEqual(a["distinct_moth_positive_plant_names"], 1)
        self.assertEqual(a["moth_plant_island_name_rows"], 2)

    def test_no_moth_record_is_not_confirmed_biological_absence(self):
        g = self.globi[self.globi.guild.eq("butterfly_family_visit")]
        a = m.scenario(g, self.flora, self.photos)
        self.assertEqual(a["distinct_moth_positive_plant_names"], 0)
        self.assertEqual(a["islands_with_at_least_1_moth_linked_plant"], 0)

    def test_quarterly_viirs_not_annual(self):
        a = m.photo_coverage(self.photos)
        self.assertIsNone(a["annual_VIIRS_eligible_island_years"])
        self.assertEqual(a["quarterly_proxy_light_islands_NOT_annual"], 1)
        self.assertEqual(a["positive_GloBI_moth_tagged_PLANT_photos_NOT_moth_visits"], 4)

    def test_source_validation_refuses_inconsistent_cohort(self):
        with self.assertRaisesRegex(ValueError, "Unexpected GloBI source frame"):
            m.validate(self.globi, self.flora, self.photos)


if __name__ == "__main__":
    unittest.main()
