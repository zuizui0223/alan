"""Synthetic fail-closed checks for source-only moth visitor/PCNL inventory."""
from __future__ import annotations

import math
import sqlite3
import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from scripts.analysis import alan_moth_species_historical_pcnl_inventory_20261109 as mod


class TestMothIslandPCNL(unittest.TestCase):
    def setUp(self):
        self.islands = gpd.GeoDataFrame(
            {"island": ["GSHHG-test-A"]},
            geometry=[box(-0.5, -0.5, 0.5, 0.5)],
            crs="EPSG:4326")

    def record(self, photo, obs, lat="0.2", lon="0.2", year="2024",
               taxon="10", db="true", walton="false", butterfly="true"):
        return {
            "photo_id": photo, "observation_uuid": obs,
            "inat_taxon_id": taxon, "observed_on": year + "-06-15",
            "latitude": lat, "longitude": lon,
            "valid_observed_date": "true", "valid_geopoint": "true",
            "moth_database_pollinates": db,
            "walton2020_nocturnal_pollen_carriage": walton,
            "butterfly_record_positive": butterfly,
        }

    def test_roster_validity_without_inventing_negative_guilds(self):
        rows = [self.record("1", "u1"), self.record("2", "u1"),
                self.record("3", "u2", year="2022"),
                self.record("4", "u3", taxon="unknown")]
        x = mod.normalized_frame(pd.DataFrame(rows), {"10": "Test plant"})
        self.assertEqual(len(x), 2)
        self.assertEqual(set(x.species), {"Test plant"})

    def test_photo_vs_unique_observation_and_butterfly_overlap(self):
        x = mod.normalized_frame(pd.DataFrame([
            self.record("1", "same"), self.record("2", "same"),
            self.record("3", "other", butterfly="false"),
            self.record("4", "outside", lat="1.5"),
        ]), {"10": "Test plant"})
        seen = set()
        rows, n, ambiguous = mod.aggregate_source_chunk(x, self.islands, seen)
        self.assertEqual(n, 3)
        self.assertEqual(ambiguous, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][5:10], (3, 2, 3, 0, 2))
        self.assertEqual(seen, {"same", "other"})
        # Already counted the observation UUIDs in previous shards.
        rows2, _, _ = mod.aggregate_source_chunk(
            x.iloc[:1].copy(), self.islands, seen)
        self.assertEqual(rows2[0][6], 0)

    def test_exact_polygon_membership_and_ambiguity_gate(self):
        islands = gpd.GeoDataFrame({"island": ["A", "B"]},
                 geometry=[box(-0.5,-0.5,0.5,0.5), box(0,0,1,1)],
                 crs="EPSG:4326")
        x = mod.normalized_frame(pd.DataFrame(
            [self.record("1", "u1"), self.record("2", "u2", lat="-0.2", lon="-0.2")]),
            {"10": "Test plant"})
        matched, ambiguous = mod.join_geography(x, islands)
        self.assertEqual(ambiguous, 1)
        self.assertEqual(matched.photo_id.tolist(), ["2"])

    def test_sql_duplicate_source_complement_is_not_extra_species(self):
        db = sqlite3.connect(":memory:")
        db.execute(mod.TABLE)
        photo = ("A",1,2,2024,"Plant1",3,2,3,0,1,0.2,0.2)
        db.execute(mod.UPSERT, photo)
        db.execute(mod.UPSERT, ("A",1,2,2024,"Plant1",1,1,1,0,0,0.2,0.2))
        db.execute(mod.UPSERT, ("A",1,2,2024,"Plant2",1,1,0,0,1,0.2,0.2))
        summary = db.execute(
           "SELECT COUNT(DISTINCT species),SUM(n_photos),SUM(n_observations) FROM taxa"
        ).fetchone()
        self.assertEqual(summary, (2,5,4))

    def test_historic_light_excludes_observation_year_and_is_not_zero_filled(self):
        import rasterio
        from rasterio.transform import from_origin
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            for year in range(1992, 2025):
                data = np.ones((2,2), dtype="float32")
                if year == 2024:
                    data[:] = 3
                if year == 2010:
                    data[0,0] = -9999
                with rasterio.open(root / f"PCNL{year}.tif", "w",
                    driver="GTiff", height=2, width=2, count=1,
                    dtype="float32", crs="EPSG:4326",
                    transform=from_origin(0,2,1,1), nodata=-9999) as w:
                    w.write(data, 1)
            frame = pd.DataFrame([
                {"year":2019,"lon":1.5,"lat":1.5},
                {"year":2025,"lon":1.5,"lat":1.5},
                {"year":2025,"lon":0.5,"lat":1.5},
            ])
            result = mod.annual_exposure(frame, root)
            self.assertAlmostEqual(result.iloc[0][
                "pcnl_1992_to_previous_year_sum_log1p"], 27*math.log(2), places=6)
            self.assertAlmostEqual(result.iloc[1][
                "pcnl_1992_to_previous_year_sum_log1p"], 32*math.log(2)+math.log(4), places=6)
            self.assertTrue(np.isnan(result.iloc[2][
                "pcnl_1992_to_previous_year_sum_log1p"]))
            self.assertEqual(result.iloc[2].pcnl_years_observed, 32)

    def test_pc_nightlight_is_required_not_substituted_by_unmatched_dmsp_viirs(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError, "Missing annual PCNL"):
                mod.annual_exposure(
                    pd.DataFrame({"year":[2025],"lon":[0.2],"lat":[0.2]}),
                    Path(t))


if __name__ == "__main__":
    unittest.main()
