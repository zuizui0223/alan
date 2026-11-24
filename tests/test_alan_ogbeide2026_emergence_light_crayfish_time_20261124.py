"""The same 16 true stream-riparian mesocosms must remain paired across weeks."""
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
MOD=ROOT/"scripts/analysis/alan_ogbeide2026_emergence_light_crayfish_time_20261124.py"
spec=importlib.util.spec_from_file_location("ogbeide_source_week",MOD)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestOgbeideSourceWeekInteraction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.arms,cls.blob,cls.week4=m.source_verify_and_read()
        cls.r=m.data_summary(cls.arms,permutation_count=1000)

    def test_original_source_identity_and_exposure(self):
        self.assertEqual(self.blob,"a774508cdbe30641e97de61a7f55d437e15b9aec")
        self.assertEqual(sum(len(v) for v in self.arms.values()),16)
        self.assertEqual({k:len(v) for k,v in self.arms.items()},{x:4 for x in m.CONDITIONS})
        for arm in self.arms.values():
            self.assertTrue(all(len(r)==3 and r[1]>0 and r[2]>0 for r in arm))
        self.assertAlmostEqual(self.week4["no_crayfish"],.6862252902416064,places=10)
        self.assertAlmostEqual(self.week4["with_crayfish"],1.417946645109135,places=10)
        self.assertTrue(self.r["both_timepoints_after_light_on"])

    def test_all_means_and_strict_repeated_flume_interaction(self):
        z=self.r
        a=z["arm_mean_emergence_week1_week4"]
        self.assertAlmostEqual(a["alan"][0],83.46875,places=10)
        self.assertAlmostEqual(a["control"][0],94.4895833333333,places=10)
        self.assertAlmostEqual(a["cray"][0],62.4895833333333,places=10)
        self.assertAlmostEqual(a["alan+cray"][0],52.53125,places=10)
        self.assertAlmostEqual(a["alan"][1],22.78125,places=10)
        self.assertAlmostEqual(a["control"][1],33.1979166666667,places=10)
        self.assertAlmostEqual(a["cray"][1],12.8854166666667,places=10)
        self.assertAlmostEqual(a["alan+cray"][1],18.2708333333333,places=10)
        self.assertAlmostEqual(z["light_by_crayfish_interaction_RRR"]["Week 1"],.9516343950785978,places=10)
        self.assertAlmostEqual(z["light_by_crayfish_interaction_RRR"]["Week 4"],2.0662990205522886,places=10)
        self.assertAlmostEqual(z["triple_interaction_ratio_of_RRRs"],2.1713160340128606,places=10)

    def test_pairwise_covariance_and_approximate_ci(self):
        z=self.r
        self.assertAlmostEqual(z["three_way_within_flume_delta_sampling_se"],.573380741991179,places=10)
        self.assertAlmostEqual(z["welch_satterthwaite_df"],11.916959699423794,places=8)
        ci=z["exploratory_paired_flume_t95_triple_ratio"]
        self.assertAlmostEqual(ci[0],.62193353,places=7)
        self.assertAlmostEqual(ci[1],7.58057425,places=7)
        self.assertAlmostEqual(z["exploratory_two_sided_t_p"],.20141110528344083,places=8)
        self.assertTrue(ci[0]<1<ci[1])
        self.assertEqual(z["independent_meta_k_increment"],0)

    def test_permutation_is_conditional_sensitivity_not_false_significance(self):
        z=self.r["conditional_2_by_2_random_label_permutation_sharp_null"]
        self.assertEqual(z["permutations_with_equal_cell_sizes"],1000)
        self.assertGreater(z["approx_two_sided_p"],.15)
        self.assertLess(z["approx_two_sided_p"],.30)
        with self.assertRaises(ValueError):m.data_summary(self.arms,permutation_count=30)

if __name__=="__main__":
    unittest.main()
