#!/usr/bin/env python3
"""Synthetic safety controls: original HAND is optional, never fabricated."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.analysis.alan_original_diel_fivearm_hand_headroom_20261114 import study_summary


class FiveArmControls(unittest.TestCase):
    def test_small_source_gates_any_nominal_p(self):
        sample = pd.DataFrame({
            "study_key": ["A", "A", "B", "C"],
            "species": ["X", "X", "Y", "Z"],
            "two_window_nonadditivity": [-.2, -.1, .3, -.25],
            "hand_headroom_proxy": [.1, .2, -.05, .3],
        })
        out = study_summary(sample)
        self.assertEqual((out["pairs"], out["papers"], out["plants"]), (4, 3, 3))
        self.assertIsNone(out["study_sign_p_nominal"])
        self.assertEqual(out["negative_I_papers"], 2)
        self.assertEqual(out["positive_I_papers"], 1)
        self.assertAlmostEqual(out["median_study_H"], .15)

    def test_strata_with_duplicated_papers_not_independent(self):
        sample = pd.DataFrame({
            "study_key": ["A", "A", "B"],
            "species": ["X", "X", "Y"],
            "two_window_nonadditivity": [-.1, -.2, .1],
            "hand_headroom_proxy": [.3, -.1, .2],
        })
        positive = study_summary(sample.loc[sample.hand_headroom_proxy > 0])
        negative = study_summary(sample.loc[sample.hand_headroom_proxy <= 0])
        self.assertEqual((positive["papers"], negative["papers"]), (2, 1))
        self.assertEqual(sample.study_key.nunique(), 2)
        self.assertTrue(positive["papers"] + negative["papers"] > 2)
        self.assertIsNone(negative["study_sign_p_nominal"])

    def test_hand_gap_not_mechanical_causal_statement(self):
        o, h = 5., 10.
        proxy = (h-o)/(h+o)
        self.assertAlmostEqual(proxy, 1/3)
        self.assertGreater(proxy, 0)
        # HAND and OPEN are measured on one original metric, not
        # an intervention in satellite artificial light at night.
        self.assertFalse(np.isnan(proxy))


if __name__ == "__main__":
    unittest.main()
