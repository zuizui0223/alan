"""Exact conditional-within-block Lu co-occurrence test."""
import math
import unittest
from scripts.analysis.alan_lu2026_sameflower_conditional_overlap_20261120 import (
  hypergeom_pmf, exact_overlap
)

class ConditionalOverlap(unittest.TestCase):
    def test_exact_hypergeometric_support_and_normalization(self):
        z=hypergeom_pmf(5,3,2)
        self.assertEqual(sorted(z),[0,1,2])
        self.assertAlmostEqual(sum(z.values()),1.)
        self.assertAlmostEqual(sum(k*v for k,v in z.items()),3*2/5.)

    def test_full_overlap_when_both_margins_cover_all(self):
        z=hypergeom_pmf(4,4,2)
        self.assertEqual(z,{2:1.})

    def test_extreme_overlap_is_fully_one_sided_under_null(self):
        rows=[]
        for b in ("A","B","C","D","E"):
            for j in range(4):
                rows.append({"block":b,"plant":f"{b}_{j}",
                  "moth_loss":j<2,
                  "day_fly_gain":j<2,
                  "night_fly_gain":j%2==0,
                  "either_fly_gain":j<2})
        z=exact_overlap(rows,"day_fly_gain")
        self.assertEqual(z["observed_moth_loss_and_fly_gain_cells"],10)
        self.assertAlmostEqual(z["null_expected_intersection"],5)
        self.assertEqual(z["independent_original_spatial_light_pairs"],5)
        self.assertAlmostEqual(z["exact_one_sided_enrichment_probability"],(1/6)**5)
        self.assertLess(z["exact_probability_ordering_two_sided_p"],.01)

    def test_independent_blocks_with_uninformative_margins(self):
        rows=[]
        for j in range(4):
            rows.append({"block":"B1","moth_loss":j<2,
                         "either_fly_gain":True})
        z=exact_overlap(rows,"either_fly_gain")
        self.assertEqual(z["observed_moth_loss_and_fly_gain_cells"],2)
        self.assertAlmostEqual(z["null_expected_intersection"],2)
        self.assertAlmostEqual(z["exact_one_sided_enrichment_probability"],1.)

if __name__=="__main__":
    unittest.main()
