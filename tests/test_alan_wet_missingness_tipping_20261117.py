"""Source-locked mathematical missingness tipping diagnostics, not a data imputation."""
import unittest
from scripts.analysis import alan_wet_moisture_missingness_tipping_20261117 as t


class TestWetMissingnessTipping(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=t.audit()

    def test_unchanged_source_and_denominators(self):
        r=self.r
        self.assertEqual(r["source_species_total"],3901)
        self.assertEqual(r["matched_original_strata"],49)
        self.assertEqual(r["known_case_count"],158)
        self.assertEqual(r["matched_strata_case_total"],257)
        self.assertEqual(r["matched_strata_case_unknown"],99)
        self.assertEqual(r["matched_strata_reference_total"],2797)
        self.assertEqual(r["matched_strata_reference_unknown"],1378)
        self.assertAlmostEqual(r["observed_complete_case_difference"],.5922104455502863)

    def test_hypothetical_shifts_zero_the_contrast(self):
        r=self.r
        d=r["same_role_missing_equals_observed_mean_scenario"]
        a=r["weight_case_missing_alpha"]
        b=r["weight_reference_missing_beta"]
        self.assertGreater(d,0)
        self.assertAlmostEqual(d-a*r["tipping_case_unknown_points_below_known_case"],0)
        self.assertAlmostEqual(d-b*r["tipping_reference_unknown_points_above_known_reference"],0)
        self.assertAlmostEqual(d-(a+b)*r["tipping_symmetric_each_side_points"],0)

    def test_no_meta_or_mechanism_promotion(self):
        self.assertIn("no missing value is actually imputed",self.r["scope"])
        self.assertIn("26 candidate species",self.r["excluded_unmatched_strata"])


if __name__=="__main__":
    unittest.main()
