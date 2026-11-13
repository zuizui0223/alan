#!/usr/bin/env python3
"""Original day/night source controls: true paired experiment site, time and unit."""
import sys, unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_DvN_fieldsite_year_source_20261113 import (
    numerical_weighted_mean, KEY, PAIRS, gitblob, STUDY_BLOB
)
class DielExperimentSourceTests(unittest.TestCase):
    def test_observed_method_grain_includes_site_and_YEAR(self):
        for k in ("study_ID","plant_species","treatment_effectiveness_metric",
                  "coordinates_lat","coordinates_lon","year_start","month_start","effectiveness_units"):
            self.assertIn(k,KEY)
        self.assertEqual(len(PAIRS),2)
    def test_weighted_pair_mean_uses_actual_sample_sizes(self):
        x=pd.DataFrame({"effectiveness_value":["0","2","4"],"sample_size":["5","2","3"]})
        self.assertAlmostEqual(numerical_weighted_mean(x),1.6)
        x.loc[0,"effectiveness_value"]="-2"
        self.assertIsNone(numerical_weighted_mean(x))
    def test_author_blob_guard_rejects_dummy_data(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as folder:
            p=Path(folder)/"test.csv"
            p.write_text("not the original author raw study")
            self.assertNotEqual(gitblob(p),STUDY_BLOB)
if __name__=="__main__":unittest.main()
