"""Exact-site wet/dry floral contrast, no matching unobserved arms to zero."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_giavi2021_same_meadow_wetdry_floral_20261120 import event_gaps

class SameMeadowSourceComparisons(unittest.TestCase):
    def test_same_event_light_minus_dark_double_difference(self):
        rows=[]
        # Same exact F values with positive floral cover, two meadows.
        for site in ("1","6"):
            for repeat in ("a","b"):
                for species,trait,dark,lit in [
                    ("wet",8.,12.,6.),
                    ("dry",4.,8.,8.),
                    ("medium",5.,5.,5.)]:
                    rows.append({"site":site,"pair_key":site+"|"+repeat,
                                 "plant":species,"F_PLANTATT":trait,
                                 "flowers_C":dark,"flowers_T":lit,
                                 "delta_flower":np.log1p(lit)-np.log1p(dark)})
        x=event_gaps(pd.DataFrame(rows),"F_PLANTATT")
        self.assertEqual(x["original_matched_same_site_same_sampling_event_comparisons"],4)
        self.assertEqual(x["independent_original_meadow_site_pairs_with_qualifying_event"],2)
        self.assertEqual(x["source_direction_negative_meadow_pairs"],2)
        self.assertAlmostEqual(x["exact_two_sided_signflip_min_resolution_if_nonzero_all_sites"],.5)
        self.assertTrue(x["cannot_infer_wet_group_specific_light_effect"])

    def test_noncoobserved_meadow_not_counted_as_biological_zero(self):
        rows=[
         dict(site="1",pair_key="1|a",plant="wet",F_PLANTATT=8.,
              flowers_T=2.,flowers_C=6.,delta_flower=np.log(3)-np.log(7)),
         dict(site="2",pair_key="2|a",plant="dry",F_PLANTATT=4.,
              flowers_T=7.,flowers_C=8.,delta_flower=np.log(8)-np.log(9))
        ]
        z=event_gaps(pd.DataFrame(rows),"F_PLANTATT")
        self.assertEqual(z["original_matched_same_site_same_sampling_event_comparisons"],0)
        self.assertEqual(z["independent_original_meadow_site_pairs_with_qualifying_event"],0)
        self.assertIsNone(z["source_within_meadow_light_minus_dark_group_difference_site_mean"])
if __name__=="__main__":
    unittest.main()
