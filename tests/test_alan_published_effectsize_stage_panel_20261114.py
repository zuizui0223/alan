"""Targeted tests; synthetic values are not additional ecological observations."""
import unittest, math
from scripts.analysis.alan_published_effectsize_stage_panel_20261114 import (
    log_rr,bounded_sd_ok,nonnegative_sd_ok,aggregate,re_meta)

class StageMetaTests(unittest.TestCase):
    def test_ratio_direction(self):
        y,v=log_rr(20,4,10,10,3,10)
        self.assertAlmostEqual(y,math.log(2));self.assertAlmostEqual(v,.013)
    def test_missing_zero_mean_not_invented(self):
        with self.assertRaises(ValueError):log_rr(0,4,10,10,3,10)
    def test_reject_invalid_n(self):
        with self.assertRaises(ValueError):log_rr(20,4,1,10,3,10)
    def test_bounded_percentage_rejects_source_sd(self):
        self.assertFalse(bounded_sd_ok(35.99,146.3282884,20,0,100))
    def test_plausible_percentage_sd(self):
        self.assertTrue(bounded_sd_ok(50,20,20,0,100))
    def test_bounded_proportion_same_rule(self):
        self.assertFalse(bounded_sd_ok(.36,1.46,20,0,1))
    def test_full_dependence_duplicate_not_extra_information(self):
        y,v=aggregate([.5,.5],[.04,.04],1)
        self.assertAlmostEqual(y,.5);self.assertAlmostEqual(v,.04)
    def test_independence_variance(self):
        y,v=aggregate([.5,.5],[.04,.04],0)
        self.assertAlmostEqual(v,.02)
    def test_bad_rho_rejected(self):
        with self.assertRaises(ValueError):aggregate([1],[1],-1)
    def test_meta_constant_and_conservative_interval(self):
        z=re_meta([.2]*5,[.04]*5)
        self.assertAlmostEqual(z['estimate'],.2)
        self.assertAlmostEqual(z['tau2'],0,places=6)
        self.assertGreater(z['ci95'][1]-.2,1.96*math.sqrt(.04/5))
    def test_nonnegative_upper_sd_bound(self):
        self.assertFalse(nonnegative_sd_ok(.15,1.055082935,23))
        self.assertFalse(nonnegative_sd_ok(.3,1.9,23))
        self.assertTrue(nonnegative_sd_ok(1,2,4))
    def test_nonnegative_all_zero(self):
        self.assertTrue(nonnegative_sd_ok(0,0,4))
        self.assertFalse(nonnegative_sd_ok(0,.1,4))
    def test_one_study_no_generalized_inference(self):
        with self.assertRaises(ValueError):re_meta([.2],[.04])

if __name__=='__main__':unittest.main()
