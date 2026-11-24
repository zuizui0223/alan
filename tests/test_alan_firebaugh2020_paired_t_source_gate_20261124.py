"""Rounding-aware Firebaugh 2020 four physical pairs audit."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"database/alan_experiments/arthropod_expansion_20261122/cross_programme/firebaugh2020_paired_t_source_gate_20261124.py"
spec=importlib.util.spec_from_file_location("firebaugh2020",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestAuthorStatistics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=m.compute()

    def test_only_four_randomized_pairs_across_both_years(self):
        x=self.result
        self.assertEqual(x["n_true_independent_randomized_plot_pairs"],4)
        self.assertTrue(x["two_annual_outcomes_are_dependently_measured"])
        self.assertEqual(x["author_reported_total_collected"],9925)
        self.assertEqual(x["independent_meta_programmes_added"],0)
        self.assertFalse(x["quantitative_meta_eligible_today"])

    def test_2015_source_paired_t_and_not_an_lnRR(self):
        x=self.result["source_rounded_paired_differences"]["2015"]
        self.assertAlmostEqual(x["implied_SE_difference_from_rounded_t3"],2.0899653979238755,places=12)
        self.assertAlmostEqual(x["approx_paired_t3_95_difference"][0],-.6112026587947916,places=12)
        self.assertAlmostEqual(x["approx_paired_t3_95_difference"][1],12.691202658794792,places=12)
        self.assertAlmostEqual(x["source_reported_lit_to_dark_relative_description"],1.61)
        self.assertTrue(x["same_four_pairs_as_other_year"])

    def test_2016_source_paired_t_and_broad_interval(self):
        x=self.result["source_rounded_paired_differences"]["2016"]
        self.assertAlmostEqual(x["implied_SE_difference_from_rounded_t3"],2.142857142857143,places=12)
        self.assertAlmostEqual(x["approx_paired_t3_95_difference"][0],-6.969527797037706,places=12)
        self.assertAlmostEqual(x["approx_paired_t3_95_difference"][1],6.669527797037706,places=12)

    def test_year_contrast_requires_original_cross_year_covariance(self):
        x=self.result
        self.assertAlmostEqual(x["difference_2016_minus_2015_in_pair_mean_differences"],-6.19)
        low,high=x["unknown_year_to_year_covariance_SE_bounds"]
        self.assertAlmostEqual(low,.052891744933267315,places=12)
        self.assertAlmostEqual(high,4.232822540781019,places=12)
        self.assertTrue(x["cannot_infer_year_by_light_interaction_p_without_joint_covariance"])
        self.assertTrue(x["author_trophic_group_interaction"]["not_lnRR"])
        self.assertEqual(x["author_trophic_group_interaction"]["model_coefficient"],.26)
        self.assertEqual(x["author_trophic_group_interaction"]["model_SE"],.46)

if __name__=="__main__":
    unittest.main()
