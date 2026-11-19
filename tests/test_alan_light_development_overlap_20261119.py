"""Synthetic tests: ALAN matching must never use botanical outcome to form pairs."""
import importlib.util
import pathlib
import unittest

import numpy as np
import pandas as pd

S = pathlib.Path(__file__).resolve().parents[1] / "scripts/analysis/alan_light_development_overlap_20261119.py"
spec = importlib.util.spec_from_file_location("overlap", S)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def make_source():
    rows = []
    for block in range(15):
        for j in range(8):
            bright = (j % 2 == 0)
            rows.append(dict(unit_id=1000*block+j, block=block,
              light=2+j*.1 if bright else -2-j*.1,
              built=.20 + .01*j, built0=2.+.1*j,
              light0=10.+.1*j, start=1999, span=22,
              temp_trend=.1, rain_trend=.1,
              effort=.02, days=.02, spring=0., summer=0.,
              freq0=10.+.1*j, area0=10.+.1*j,
              freq_trend=.1, area_trend=.1,
              candidate_rate=-1. if bright else 0.,
              reference_rate=0., rate_difference=-1. if bright else 0.))
    return pd.DataFrame(rows)

class TestMatching(unittest.TestCase):
    def test_pairs_ignore_botanical_outcome(self):
        a = make_source()
        p, f, s = m.pair_sites(a)
        self.assertGreaterEqual(len(p), 30)
        self.assertEqual(p.block.nunique(), 15)
        b = a.copy()
        b["rate_difference"] = np.arange(len(b), dtype=float)[::-1]
        b["candidate_rate"] = b.rate_difference
        q, _, _ = m.pair_sites(b)
        for column in ("block", "light_difference", "d_built", "d_light0"):
            np.testing.assert_allclose(p[column].to_numpy(), q[column].to_numpy())

    def test_block_and_built_caliper_stop(self):
        a = make_source()
        a.loc[a.light.gt(0.), "built"] += 100.
        p, _, _ = m.pair_sites(a)
        self.assertEqual(len(p), 0)

    def test_identical_pairs_nested_water_regression(self):
        a = make_source()
        p, _, _ = m.pair_sites(a)
        z = m.same_cohort_water_adjustment(p)
        self.assertEqual(z["same_matched_pairs"], len(p))
        self.assertTrue(z["both_models_include_baseline_water_and_climate"])
        self.assertEqual(z["status"], "POSTHOC_SENSITIVITY_ONLY")

    def test_report_has_no_causal_claim(self):
        a = make_source()
        p, f, s = m.pair_sites(a, water_matched=True)
        summary = m.summarize_pairs(p, f, s)
        self.assertEqual(summary["matched_pairs"], len(p))
        self.assertIn("postselected", summary["interpretation"])
        self.assertIn(summary["status"], ("SOURCE_SUPPORT_STOP",
              "EXPLORATORY_COMMON_SUPPORT"))

if __name__ == "__main__":
    unittest.main()
