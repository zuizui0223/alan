import unittest
from scripts.analysis import alan_A_Boom2020_original_SE_covariance_20261115 as m
class TestBoomCovariance(unittest.TestCase):
    def test_original_source_constants(self):
        self.assertEqual(m.EXPECTED_SE,{"Dark":40.0,"White":28.0,"Green":20.0,"Red":28.0})
        self.assertEqual(m.EXPECTED_MEAN["Dark"],117)
        self.assertEqual(m.FIELD_SITES,4)
    def test_correlation_bound(self):
        x=m.covariance_bounds(117,40,61,20)
        self.assertEqual(x["reported_mean_treatment_minus_dark"],-56)
        self.assertEqual(x["possible_pair_difference_SE_min_if_rho_plus1"],20)
        self.assertEqual(x["possible_pair_difference_SE_max_if_rho_minus1"],60)
        self.assertAlmostEqual(x["hypothetical_best_case_abs_t_if_four_independent_matched_sites"],2.8)
        self.assertFalse(x["hypothetical_best_case_t_exceeds_df3_threshold"])
        self.assertAlmostEqual(x["illustrative_independent_arm_SE_if_rho_zero"],(40**2+20**2)**.5)
    def test_unknown_covariance_not_meta(self):
        x=m.covariance_bounds(117,40,86,28)
        self.assertTrue(x["published_covariance_rho_not_available"])
        self.assertTrue(x["not_a_computed_meta_lnRR_sampling_variance"])
        self.assertEqual(x["original_full_4_treatment_ANOVA_p"],.63)
    def test_reject_nonpositive_seed_mean(self):
        with self.assertRaises(ValueError):m.covariance_bounds(0,40,61,20)
    def test_multiple_spectra_one_source(self):
        self.assertEqual(len(m.SOURCE_NAMES),4)
if __name__=="__main__":unittest.main()
