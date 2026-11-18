"""Guard original mesocosm randomization and avoid pseudo-replication."""
import itertools
import math
import unittest
import numpy as np
from scripts.analysis.alan_true_parkinson2020_aquatic_spider_20261118 import (
    exact_unblocked_test,bootstrap_mean_ratio,bh,PRESET_PRIMARY,SOURCE_BLOB_SHA1
)

class TestRealExperimentInferenceUnit(unittest.TestCase):
    def test_original_20_randomization_labels_not_60_independent_dates(self):
        treated=np.array([True]*10+[False]*10)
        values=np.array([2.0]*10+[1.0]*10)
        out=exact_unblocked_test(values,treated)
        self.assertEqual(out["enumerated_unblocked_label_sets"],math.comb(20,10))
        self.assertAlmostEqual(out["observed_raw_mean_difference"],1.0)
        self.assertAlmostEqual(out["descriptive_unblocked_fixed_10_per_arm_permutation_p"],
                               2/math.comb(20,10))
        self.assertTrue(out["not_exact_for_original_unobserved_fish_blocking"])

    def test_unchanged_source_record_fields_and_four_endpoint_gate(self):
        self.assertEqual(len(PRESET_PRIMARY),4)
        self.assertEqual(set(SOURCE_BLOB_SHA1),{
            "Spider and Prey Insect Abundance.csv",
            "Aquatic and Terrestrial Invertebrate Inputs.csv",
            "Spider Body Mass.csv"
        })
        self.assertEqual(len(set(SOURCE_BLOB_SHA1.values())),3)

    def test_bootstrap_resamples_same_mesocosm_not_body_mass(self):
        rng=np.random.default_rng(20261118)
        out=bootstrap_mean_ratio([4.]*10,[2.]*10,rng,n=500)
        self.assertEqual(out["bootstrap_n"],500)
        self.assertEqual(out["effective_bootstrap_comparisons"],500)
        self.assertAlmostEqual(out["percentile_95CI"][0],2.)
        self.assertAlmostEqual(out["percentile_95CI"][1],2.)

    def test_bh_is_not_unadjusted_signal_fishing(self):
        q=bh([0.001,0.02,0.04,0.8])
        self.assertTrue(all(0<=x<=1 for x in q))
        self.assertAlmostEqual(q[0],.004)
        self.assertAlmostEqual(q[-1],.8)
        self.assertAlmostEqual(q[2],.04*4/3)

if __name__=="__main__":
    unittest.main()
