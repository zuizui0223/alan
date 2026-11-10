import unittest
import pandas as pd
from scripts.acquisition.alan_tier09_vs_full112428_measured_geography_20261110 import (
    require, boolean, evaluate
)

class PrecolourFullSourceJoinTests(unittest.TestCase):
    def test_fail_closed_binary_label(self):
        with self.assertRaises(ValueError):
            boolean(pd.Series(["true","unknown"]),"roi")

    def test_exact_photo_id_overlap_rejected_before_modelling(self):
        # Make a minimal source-ID table; fail before the fixed historical
        # cohort checks, so a future pipeline cannot silently double-measure.
        photos=pd.DataFrame({
            "photo_id":["shared"],"real_numeric_ROI9":[True],
            "target_id":["2024|a|i1"],"year":["2024"],"land_id":["i1"]})
        new=pd.DataFrame({"photo_id":["shared"]})
        with self.assertRaises(ValueError):
            evaluate(photos,pd.DataFrame(),new,pd.DataFrame(),{})

    def test_source_keys_normalized_as_strings(self):
        x=pd.DataFrame({"land_id":["0004"],"year":["2024"],"cell_id":["1:2"]})
        self.assertEqual(
            list(x[["land_id","year","cell_id"]].astype(str)
                 .itertuples(index=False,name=None)),
            [("0004","2024","1:2")])

if __name__=="__main__":
    unittest.main()
