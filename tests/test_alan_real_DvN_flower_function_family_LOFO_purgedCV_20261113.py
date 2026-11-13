#!/usr/bin/env python3
"""Ensure independent family holdout cannot leak same plant to training."""
import sys
import unittest
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_real_DvN_flower_function_family_LOFO_purgedCV_20261113 import (
    lofo_models, purged_study_out, N_PAIRS,N_STUDIES,N_SPECIES,N_FAMILIES, ALPHA
)

class DielFunctionOriginalFamilyTests(unittest.TestCase):
    def test_frozen_original_trait_cohort(self):
        self.assertEqual((N_PAIRS,N_STUDIES,N_SPECIES,N_FAMILIES),(196,73,74,28))
        self.assertEqual(ALPHA,5.)

    def test_no_fake_out_of_family_CV_without_many_sources(self):
        x=pd.DataFrame({"study_key":["a","a","b","b"],"plant_family":["X","Y","X","Y"],
                        "genus":["P","Q","P","Q"],"night_bloom":[1,0,1,0]})
        with self.assertRaises(ValueError):
            lofo_models(x,"plant_family")

    def test_incomplete_study_leakage_group_refuses_fake_CV(self):
        # Missing 196 authenticated source outcomes is a hard error.
        df=pd.DataFrame({
            "study_key":["a","b"],
            "species":["x","y"],"night_bloom":[1,0],
            "plant_family":["Fam1","Fam2"]})
        with self.assertRaises(ValueError):
            purged_study_out(df)

if __name__=="__main__":
    unittest.main()
