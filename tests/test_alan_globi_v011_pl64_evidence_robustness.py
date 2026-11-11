"""Synthetic mathematical safeguards, never additional biological observations."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts" / "analysis"))
import alan_globi_v011_pl64_evidence_robustness as m


class PL64SourceMath(unittest.TestCase):
    def test_original_equal_species_moderator_algebra(self):
        names = ["Plant a", "Plant b", "Plant c"]
        samples = pd.DataFrame({"Plant a": [1, 0, 1], "Plant b": [0, 1, 1],
                                "Plant c": [1, 0, 0]})
        PL = pd.Series([0., 1., 2.], index=names)
        result = m.weighted_species_response(samples, np.array([3, 4, 5]), PL)
        self.assertTrue(np.allclose(result, [-100/3, 0, -100/5]))
        self.assertAlmostEqual(float(result.sum()), -100/3 - 100/5)

    def test_direct_moderator_contrast_linear(self):
        cols = ["a", "b", "c", "d", "e", "f"]
        data = pd.DataFrame(np.array([[1, 0, 1, 0, 0, 1],
                        [0, 1, 0, 1, 0, 0], [1, 1, 1, 1, 1, 1]]), columns=cols)
        pl1 = pd.Series([0., 1., 2.], index=cols[:3])
        pl2 = pd.Series([0., 2., 4.], index=cols[3:])
        d = np.array([3, 3, 6])
        p = m.weighted_species_response(data, d, pl1)
        q = m.weighted_species_response(data, d, pl2)
        self.assertTrue(np.allclose(p - q,
            m.weighted_species_response(data, d, pl1)
            - m.weighted_species_response(data, d, pl2)))

    def test_wls_FWL_original_population_and_cluster_blocks(self):
        # Exactly 195 two-class year×region groups and 78 repeating geography
        # blocks, matching the source design without copying real photo cells.
        rng = np.random.default_rng(971)
        rows = []
        for block in range(78):
            for yy in range(3 if block < 39 else 2):
                ry = f"b{block}|{yy}"
                idx = len({x["region_year"] for x in rows})  # synthetic only
                n_island = 5 if idx < 153 else 4
                n_mainland = 51 if idx < 15 else 50
                for land, n in (("island", n_island),
                                ("mainland", n_mainland)):
                    for j in range(n):
                        light = rng.uniform(0, 5)
                        rows.append({"region_year": ry, "region5": f"b{block}",
                          "land_class":land, "island":float(land == "island"),
                          "loglight":light, "light_island":light*float(land=="island"),
                          "logpop": rng.uniform(0, 4),
                          "logphotos": rng.uniform(1, 5),
                          "logspecies": rng.uniform(.7, 3),
                          "bio5":rng.uniform(5,30),
                          "lograin":rng.uniform(2,7),
                          "elevation_km":rng.uniform(0,2)})
        df = pd.DataFrame(rows)
        self.assertEqual(len(df), 10698)
        self.assertEqual(int(df.island.sum()), 933)
        self.assertEqual(df.region5.nunique(), 78)
        self.assertEqual(df.region_year.nunique(), 195)
        y = rng.normal(0, 1, len(df)) + 0.35*df.light_island
        fit = m.make_fit(df)
        beta = fit(y.to_numpy())
        grp = df.groupby(["region_year","land_class"]).size()
        index = pd.MultiIndex.from_frame(df[["region_year","land_class"]])
        w = grp.groupby(level=0).min().reindex(df.region_year).to_numpy(float)/grp.reindex(index).to_numpy(float)
        df["outcome"] = y
        expected = smf.wls("outcome ~ island + loglight + light_island + logpop + logphotos + logspecies + bio5 + lograin + elevation_km + C(region_year)",
                           data=df, weights=w).fit()
        self.assertAlmostEqual(beta["beta"], expected.params["light_island"], places=7)


if __name__ == "__main__":
    unittest.main()
