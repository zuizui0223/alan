"""Strict original 5 spatial-block sign-flip and author coding guards."""
import unittest
from scripts.analysis.alan_lu2026_original_fiveblock_guild_replacement import (
  exact_five_signflip,normalize_time,normalize_light,need
)

class TestLuOriginalFiveBlocks(unittest.TestCase):
    def test_night_day_and_original_treatments(self):
        self.assertEqual(normalize_time("Nocturnal"),"night")
        self.assertEqual(normalize_time("Diurnal"),"day")
        self.assertEqual(normalize_light("Lit"),"lit")
        self.assertEqual(normalize_light("Unlit"),"dark")
        self.assertTrue(normalize_time("unknown").startswith("unknown:"))

    def test_five_blocks_minimum_two_sided_p_is_sixteenth(self):
        r=exact_five_signflip([-2.,-2.,-1.,-1.,-1.])
        self.assertEqual(r["number_of_source_randomization_signflips"],32)
        self.assertAlmostEqual(r["descriptive_exact_32_paired_label_signflips_two_sided_p"],.0625)
        self.assertTrue(r["not_a_large_n_significant_test"])
        self.assertEqual(r["independent_block_negative_effects"],5)

    def test_never_promote_single_plant_or_visitation_event(self):
        with self.assertRaises(ValueError):
            exact_five_signflip([-1.]*4)
        with self.assertRaises(ValueError):
            exact_five_signflip([-1.]*6)

if __name__=="__main__":
    unittest.main()
