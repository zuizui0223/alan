"""Tests for static visit replacement arithmetic and non-causal source guards."""
import unittest
from scripts.analysis.alan_page2021_replacement_threshold_20261119 import (
    calculate, bootstrap, summary
)

class VisitReplacementMath(unittest.TestCase):
    def test_same_visits_not_always_same_pollen(self):
        a=calculate([8,2],[5,1])
        self.assertAlmostEqual(a["required_new_visits_per_lost_visit"],5.)
        self.assertAlmostEqual(a["equal_visits_restored_pollen_ratio"],10./42.)
        self.assertTrue(a["requires_more_than_equal_visits"])
        self.assertAlmostEqual(a["conditional_half_top_visits_shift_pollen_ratio"],(1+10./42.)/2)

    def test_equal_efficiency_full_compensation(self):
        a=calculate([8,2],[2,2])
        self.assertAlmostEqual(a["equal_visits_restored_pollen_ratio"],1.)
        self.assertAlmostEqual(a["required_new_visits_per_lost_visit"],1.)
        self.assertFalse(a["requires_more_than_equal_visits"])

    def test_inefficient_dominant_replace_gain(self):
        a=calculate([8,2],[1,5])
        self.assertLess(a["required_new_visits_per_lost_visit"],1.)
        self.assertGreater(a["equal_visits_restored_pollen_ratio"],1.)
        self.assertGreater(a["equal_visits_functional_gain"],0.)
        self.assertEqual(a["equal_visits_functional_shortfall"],0.)

    def test_zero_alternative_efficiency_is_not_fabricated(self):
        a=calculate([8,2],[5,0])
        self.assertTrue(a["requires_infinite_or_unavailable_visits"])
        self.assertIsNone(a["required_new_visits_per_lost_visit"])
        self.assertAlmostEqual(a["equal_visits_restored_pollen_ratio"],0.)

    def test_two_guilds_not_required_to_infer_conditional_result(self):
        a=calculate([8,2,1],[2,4,8])
        self.assertEqual(a["taxa"],3)
        self.assertTrue(a["new_visits_not_actually_observed"])

    def test_cluster_sample_size_gate(self):
        a=calculate([8,2],[1,3])
        small=[dict(a,source_study="PaperA",source_plant="PlantA")]
        self.assertEqual(bootstrap(small,3)["status"],
                         "UNDER_FIVE_INDEPENDENT_STUDIES_STOP")
        self.assertEqual(summary(small)["independent_original_studies"],1)

if __name__ == "__main__":
    unittest.main()
