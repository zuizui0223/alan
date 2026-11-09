"""Synthetic guards for tier09 SHA-gated existing-cell VIIRS support only."""
import unittest
import pandas as pd
from scripts.acquisition.alan_tier09_existingVIIRS_cell_support_gate_20261110 import (
    score_support, cohort_support, LIGHT_MIN_LOG_SPAN,
)

class ExistingViirsSupportTest(unittest.TestCase):
    def fixture(self):
        rows = []
        for cell in ("1_1", "2_1"):
            for j in range(6):
                rows.append({
                    "photo_id": cell + str(j), "land_id": "45", "year": "2024",
                    "cell_id": cell, "target_id": "2024|" + cell + "|45",
                    "species": ("sp1", "sp2", "sp3")[j % 3],
                    "area_km2": "12", "moth_database_pollinates": "true",
                    "walton2020_nocturnal_pollen_carriage": "false",
                    "butterfly_record_positive": "true",
                })
        return pd.DataFrame(rows)

    def test_two_valid_original_cells_with_real_light_difference(self):
        d = self.fixture()
        v = pd.DataFrame({
            "target_id": ["2024|1_1|45", "2024|2_1|45"],
            "viirs_median": [0.0, 2.0],
            "viirs_valid_months": [4, 4],
        })
        scored = score_support(d, v)
        result, pairs = cohort_support(scored)
        self.assertEqual(result["supported_small_islands_le20km2"], 1)
        self.assertEqual(result["independent_supported_islands"], 1)
        self.assertEqual(result["qualified_existing_valid_VIIRS_sites"], 2)
        self.assertGreater(pairs.span.iloc[0], LIGHT_MIN_LOG_SPAN)

    def test_missing_grid_is_unknown_not_zero_light(self):
        d = self.fixture()
        v = pd.DataFrame({
            "target_id": ["2024|1_1|45"],
            "viirs_median": [0.0], "viirs_valid_months": [4],
        })
        scored = score_support(d, v)
        self.assertFalse(scored.loc[scored.cell_id.eq("2_1"), "existingcell_source_hit"].any())
        self.assertTrue(scored.loc[scored.cell_id.eq("2_1"), "log_light"].isna().all())
        result, _ = cohort_support(scored)
        self.assertEqual(result["independent_supported_islands"], 0)
        self.assertEqual(result["existing_VIIRS_absent_or_ineligible_qualified_sites"], 1)

    def test_identical_light_not_supported(self):
        d = self.fixture()
        v = pd.DataFrame({
            "target_id": ["2024|1_1|45", "2024|2_1|45"],
            "viirs_median": [0.5, 0.5],
            "viirs_valid_months": [4, 4],
        })
        result, _ = cohort_support(score_support(d, v))
        self.assertEqual(result["independent_supported_islands"], 0)

    def test_missing_months_not_treated_as_valid(self):
        d = self.fixture()
        v = pd.DataFrame({
            "target_id": ["2024|1_1|45", "2024|2_1|45"],
            "viirs_median": [0.0, 1.0],
            "viirs_valid_months": [4, 2],
        })
        result, _ = cohort_support(score_support(d, v))
        self.assertEqual(result["independent_supported_islands"], 0)

if __name__ == "__main__":
    unittest.main()
