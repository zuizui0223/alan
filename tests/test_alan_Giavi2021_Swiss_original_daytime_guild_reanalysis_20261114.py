#!/usr/bin/env python3
"""Swiss Giavi true paired site sign and independent-study safety tests."""
import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114 import (
    full_six_block_signflip, SOURCE_TOTAL, GUILDS, require
)
class SwissSiteLevel(unittest.TestCase):
 def test_six_site_signflip_not_1884_pseudo_replicates(self):
    x=full_six_block_signflip(np.array([-1,-2,-3,-4,-5,-6.]))
    self.assertEqual(x["site_pairs"],6)
    self.assertEqual(x["positive_site_differences"],0)
    self.assertEqual(x["negative_site_differences"],6)
    self.assertEqual(x["two_sided_full_64_site_signflip_fraction_descriptive"],.03125)
    self.assertTrue(x["not_claimed_true_assignment_randomization_test"])
 def test_original_cohort_scope(self):
    self.assertEqual(SOURCE_TOTAL,1884)
    self.assertEqual(set(GUILDS),{"Diptera","Hymenoptera","Coleoptera"})
 def test_incomplete_swiss_site_support_must_fail(self):
    with self.assertRaises(ValueError):
        full_six_block_signflip([-1,-1,1])
    with self.assertRaises(ValueError):
        full_six_block_signflip([1,2,3,4,5,float("nan")])

if __name__=="__main__":
 unittest.main()
