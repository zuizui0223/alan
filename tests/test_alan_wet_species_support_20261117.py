"""Fail-closed checks for source-based 3,901 species MOIST sensitivity."""
import math
import unittest
from scripts.analysis import alan_wet_species_support_20261117 as a


class WetSpeciesSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = a.run()

    def test_frozen_primary_not_mutated(self):
        self.assertAlmostEqual(self.r["frozen_primary_moisture_difference"],
                               0.5922104455502863)
        c = self.r["coverage"]
        self.assertEqual(c["all_species"], 3901)
        self.assertEqual(c["candidate_count"], 283)
        self.assertEqual(c["candidate_indicator_known"], 160)
        self.assertEqual(c["matched_strata_count"], 49)
        self.assertEqual(c["matched_strata_case_known"], 158)

    def test_cover_and_bins_never_recode_missing_as_zero(self):
        c = self.r["coverage"]
        self.assertLess(c["candidate_coverage"], 1.)
        self.assertLess(c["reference_coverage"], 1.)
        bins = self.r["moisture_bins_and_odds"]["bins"]
        self.assertAlmostEqual(sum(z["candidate_n"] for z in bins.values()), 158)
        self.assertAlmostEqual(sum(z["matched_ref_expected_n"] for z in bins.values()), 158)
        self.assertGreater(self.r["moisture_bins_and_odds"]["wet_gt7_conditional_MH_odds_ratio"],0)

    def test_geographic_holdout_unselected(self):
        g = self.r["geographic_sensitivity"]
        self.assertEqual(g["n_regions"], 12)
        self.assertEqual(g["positive"] + g["negative"], 12)
        self.assertEqual(g["largest_two_case_count"],94)
        self.assertEqual(g["excluding_largest_two_cases"],64)
        self.assertTrue(math.isfinite(g["excluding_largest_two_difference"]))
        self.assertEqual(len(g["leave_one"]),12)
        self.assertEqual(len(self.r["not_estimated"]),3)


if __name__=="__main__":
    unittest.main()
