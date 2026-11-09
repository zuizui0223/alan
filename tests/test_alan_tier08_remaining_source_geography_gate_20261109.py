from __future__ import annotations
import unittest
import pandas as pd

from scripts.acquisition.alan_tier08_remaining_source_geography_gate_20261109 import (
    qualified_same_year_spatial_support, summarize)


def record(island, year, site, i, light):
    return {"photo_id": f"{island}-{year}-{site}-{i}",
            "observation_id": f"uuid-{island}-{year}-{site}-{i}",
            "land_id": island, "year": str(year), "target_id": site,
            "log_light": str(light), "area_km2": "35.0",
            "wcvp_accepted_plant_name_id": str(i % 3)}


class RemainderSourceOnlyTests(unittest.TestCase):
    def test_two_sites_must_cooccur_same_year_with_true_exposure_span(self):
        rows = [record("A", 2024, "A-dark", i, 0.01) for i in range(5)]
        rows += [record("A", 2024, "A-lit", i, 0.30) for i in range(5)]
        rows += [record("B", 2023, "B-dark", i, 0.01) for i in range(5)]
        rows += [record("B", 2024, "B-lit", i, 0.30) for i in range(5)]
        result = qualified_same_year_spatial_support(pd.DataFrame(rows))
        self.assertEqual(result["source_sites_with_5_photos_3_taxa"], 4)
        self.assertEqual(result["islands_with_two_sites_sameyear_light_gap_0p05"], 1)
        self.assertEqual(result["paired_eligible_source_sites"], 2)

    def test_repeat_photos_cannot_promote_site_or_pair(self):
        rows = [record("A", 2024, "dark", i, 0.0) for i in range(5)]
        rows += [record("A", 2024, "lit", i, 0.02) for i in range(5)]
        result = qualified_same_year_spatial_support(pd.DataFrame(rows))
        self.assertEqual(result["islands_with_two_sites_sameyear_light_gap_0p05"], 0)
        rows.append(rows[0].copy())
        with self.assertRaisesRegex(ValueError, "duplicate"):
            qualified_same_year_spatial_support(pd.DataFrame(rows))

    def test_remainder_never_reuses_measured_photo_id(self):
        rows = [record("C", 2025, "s1", i, 0.2) for i in range(8)]
        novel = pd.DataFrame(rows)
        measured = novel.iloc[:3].copy()
        rest, stats = summarize(novel, measured)
        self.assertEqual(len(rest), 5)
        self.assertFalse(rest.photo_id.isin(measured.photo_id).any())
        self.assertEqual(stats["tier08_unattempted_exact_photo_ID_candidates"], 5)


if __name__ == "__main__":
    unittest.main()
