"""Synthetic precolour source invariants only; empirical counts require signed Actions artifacts."""
import unittest

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

from scripts.acquisition.alan_tier09_moth_precolour_GSHHG_photoID_gate_20261110 import (
    unique_precolour_candidates, check,
)


def fixture():
    raw = pd.DataFrame({
        "photo_id": ["11", "12", "13", "14", "15"],
        "observation_uuid": ["u1", "u1", "u2", "u3", "u4"],
        "inat_taxon_id": ["1", "1", "2", "3", "1"],
        "observed_on": ["2024-06-01", "2024-06-01", "2024-06-02",
                        "2024-06-03", "2022-06-01"],
        "valid_observed_date": ["true"] * 5,
        "valid_geopoint": ["true"] * 5,
        "frozen_loader_supported": ["true"] * 5,
        "latitude": ["1.0", "1.0", "1.1", "1.2", "1.1"],
        "longitude": ["1.0", "1.0", "1.1", "1.2", "1.1"],
        "tier": ["9"] * 5,
        "image_extension": ["jpg"] * 5,
        "photo_license": ["cc-by"] * 5,
        "moth_database_pollinates": ["true", "true", "false", "false", "false"],
        "walton2020_nocturnal_pollen_carriage": ["false"] * 5,
        "butterfly_record_positive": ["true"] * 5,
    })
    islands = gpd.GeoDataFrame(
        {"land_id": ["i1"], "area_km2": [12.0]},
        geometry=[box(0.5, 0.5, 2, 2)], crs="EPSG:4326")
    return raw, islands


class Tier09PrecolourTest(unittest.TestCase):
    def test_dedup_frozen_photos_observations_and_year(self):
        raw, islands = fixture()
        data, audit = unique_precolour_candidates(
            raw, {"13"}, {"1": "W1", "2": "W2", "3": "W3"}, islands)
        self.assertEqual(data.photo_id.tolist(), ["11", "14"])
        self.assertTrue(data.observation_uuid.is_unique)
        self.assertEqual(data.land_id.nunique(), 1)
        self.assertEqual(set(data.year), {"2024"})
        self.assertEqual(float(data.area_km2.min()), 12.0)
        self.assertEqual(audit["ambiguous"], 0)

    def test_reject_duplicate_original_photo(self):
        raw, islands = fixture()
        bad = pd.concat([raw, raw.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "Repeated"):
            unique_precolour_candidates(bad, set(), {"1": "W1"}, islands)

    def test_reject_mixed_original_tier(self):
        raw, islands = fixture()
        raw.loc[0, "tier"] = "8"
        with self.assertRaisesRegex(ValueError, "Non-tier09"):
            unique_precolour_candidates(raw, set(), {"1": "W1"}, islands)

    def test_reject_inconsistent_source_guild_flag(self):
        raw, islands = fixture()
        raw.loc[0, "moth_database_pollinates"] = "unknown"
        with self.assertRaisesRegex(ValueError, "Bad source flag"):
            unique_precolour_candidates(raw, set(), {"1": "W1", "2": "W2"}, islands)


if __name__ == "__main__":
    unittest.main()
