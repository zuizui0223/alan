"""Synthetic algebraic guards; ecological coefficients require signed 93m inputs."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_real_genus_prevalence_bias_20261110 as audit


class GenusAuditTests(unittest.TestCase):
    def test_within_island_year_sparse_effect_recovers_negative_linear_relation(self):
        y, light, obs, w, isl, iy, genus = [], [], [], [], [], [], []
        for i in range(36):
            for year in (2019, 2023, 2025):
                for j in range(4):
                    L = .1+.1*i+.3*j+.02*(year-2019)
                    genus_label = "G" + str(j % 2)
                    # Controlled genus and island×year nuisance pattern
                    y.append(.2-.07*L+.04*(j % 2)+.001*i+.002*(year-2019))
                    light.append(L)
                    obs.append(np.log1p(10+j))
                    w.append(1.)
                    isl.append(str(i))
                    iy.append(str(i)+"|"+str(year))
                    genus.append(genus_label)
        result = audit.fixed_effects_regression(
            y, light, obs, w, np.asarray(isl), [np.asarray(iy), np.asarray(genus)])
        self.assertAlmostEqual(result["beta_pp"], -7.0, places=4)
        self.assertEqual(result["n_islands"], 36)

    def test_source_missingness_and_frequency_types_are_not_moth_dependence(self):
        self.assertEqual(len(audit.GLOBI_SHA), 64)
        self.assertEqual(audit.KEY, ["island", "latbin", "lonbin", "year"])
        self.assertIn("real93m_genus", audit.SCHEMA.lower())

    def test_source_sha_rejects_noncanonical_input(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as folder:
            f = Path(folder)/"bad.csv"
            f.write_text("island,year\n1,2025\n")
            with self.assertRaisesRegex(ValueError, "Original signed"):
                audit.load_sites(f)


if __name__ == "__main__":
    unittest.main()
