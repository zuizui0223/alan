#!/usr/bin/env python3
"""Synthetic numerical tests only; live source is verified by dedicated GitHub CI."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
from alan_globi_crossland_light_robust_20261111 import fit
from alan_globi_v011_old_new_species_attribution_20261112 import (
    BETA, light_island_influence, supported_fits, classify_species)


def fake_sites():
    rng = np.random.default_rng(20261112)
    records = []
    for region in range(12):
        for year in (2023, 2024):
            for land_class in ("island", "mainland"):
                for site in range(5):
                    island = float(land_class == "island")
                    loglight = float(rng.uniform(0, 4) + island * 0.15)
                    records.append({
                        "target_id": f"t_{region}_{year}_{land_class}_{site}",
                        "key": f"k_{region}_{year}_{land_class}_{site}",
                        "region_year": f"{year}|{region}",
                        "region5": str(region), "year": year,
                        "land_class": land_class, "island": island,
                        "loglight": loglight,
                        "light_island": loglight * island,
                        "logpop": float(rng.uniform(0, 5)),
                        "logphotos": float(rng.uniform(1, 4)),
                        "logspecies": np.log1p(6),
                        "bio5": float(rng.normal(20, 3)),
                        "lograin": float(rng.uniform(0, 5)),
                        "elevation_km": float(rng.uniform(-0.1, 2)),
                        "n_species": 6
                    })
    return pd.DataFrame(records)


class TestOldNewAttribution(unittest.TestCase):
    def test_weighted_FWL_matches_existing_source_fit(self):
        sites = fake_sites()
        rng = np.random.default_rng(88)
        influence = light_island_influence(sites)
        for _ in range(4):
            counts = rng.integers(-3, 4, len(sites))
            observed = float((100 * influence * counts / sites.n_species).sum())
            original = fit(sites.assign(all_moth=counts), "all", "moth", balance=True)
            self.assertTrue(original["estimable"])
            self.assertAlmostEqual(observed, original[BETA], places=9)

    def test_additive_changed_labels_with_same_design(self):
        s = fake_sites()
        rows = []
        rng = np.random.default_rng(4)
        old_count = []
        for key in s.key:
            cl = rng.choice(["common", "old_only", "new_only", "neither"],
                            size=6, p=[0.23, 0.32, 0.10, 0.35])
            old_count.append(int(np.isin(cl, ["common", "old_only"]).sum()))
            rows.extend([{"key": key, "class": a} for a in cl])
        s["all_moth"] = old_count
        obs = pd.DataFrame(rows)
        _, responses, models = supported_fits(s, obs)
        b = lambda x: models[x][BETA]
        self.assertAlmostEqual(b("published_old"), b("old_only") + b("common"), 9)
        self.assertAlmostEqual(b("updated_v011"), b("new_only") + b("common"), 9)
        self.assertAlmostEqual(b("updated_minus_published"), b("new_only") - b("old_only"), 9)
        self.assertTrue((responses["updated_minus_published"] ==
                        responses["new_only"] - responses["old_only"]).all())

    def test_flags_preserve_taxon_status_and_do_not_impute_absence(self):
        sp = pd.DataFrame({
            "key": ["k"]*4,
            "wcvp_accepted_plant_name_id": ["1", "2", "3", "4"],
            "wcvp_accepted_name": ["Old sp", "New sp", "Both sp", "Unknown sp"]
        })
        old = pd.DataFrame({
            "guild": ["moth_family_visit", "moth_family_visit"],
            "plant_name": ["Old sp", "Both sp"]
        })
        with tempfile.TemporaryDirectory() as td:
            fp = Path(td) / "labels.csv"
            pd.DataFrame({
                "accepted_id": ["1", "2", "3"],
                "accepted_name": ["Old sp", "New sp", "Both sp"],
                "legacy_wcvp": [True, False, True],
                "v011_wcvp": [False, True, True]
            }).to_csv(fp, index=False)
            z = classify_species(sp, old, fp)
        self.assertEqual(z["class"].tolist(),
                         ["old_only", "new_only", "common", "neither"])
        self.assertFalse(z.loc[3, "new_positive"])


if __name__ == "__main__":
    unittest.main()
