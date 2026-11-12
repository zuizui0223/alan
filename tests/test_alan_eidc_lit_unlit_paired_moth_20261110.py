"""Synthetic checks of paired site inference; no biological test observations fabricated."""
from __future__ import annotations
import unittest
import pandas as pd
from scripts.analysis.alan_eidc_actual_lit_unlit_moth_pollen_20261110 import paired_site_effect

class PairInferenceTests(unittest.TestCase):
    def test_repeated_nights_do_not_increase_field_site_count(self):
        frame=pd.DataFrame([
            {"Site":11,"Lit":10,"Unlit":5},
            {"Site":11,"Lit":20,"Unlit":10},
            {"Site":21,"Lit":9,"Unlit":3},
            {"Site":21,"Lit":8,"Unlit":4},
            {"Site":30,"Lit":5,"Unlit":1}
        ])
        r=paired_site_effect(frame)
        self.assertEqual(r["independent_field_sites"],3)
        self.assertEqual(r["original_matched_sampling_events"],5)
        self.assertGreater(r["per_site_equal_exp_log_ratio"],1.)
        self.assertEqual(r["lit_more_events"],5)
    def test_distinguishes_negative_effect_and_preserves_ties(self):
        frame=pd.DataFrame([
            {"Site":11,"Lit":5,"Unlit":9},
            {"Site":21,"Lit":10,"Unlit":10},
            {"Site":30,"Lit":2,"Unlit":8},
            {"Site":40,"Lit":3,"Unlit":4}
        ])
        r=paired_site_effect(frame)
        self.assertLess(r["per_site_equal_exp_log_ratio"],1.)
        self.assertEqual(r["exact_ties"],1)
    def test_fail_closed_on_nonpairs(self):
        with self.assertRaisesRegex(ValueError,"No actual pair"):
            paired_site_effect(pd.DataFrame([{"Site":1,"Lit":4,"Unlit":3}]))
if __name__=="__main__":unittest.main()
