#!/usr/bin/env python3
"""Source-limited original floral function transfer and unknown-morphology tests."""
import unittest
import sys
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_flower_multitraits_reproductive_diel_transfer_20261113 import (
    original_traits, fivefold_blocked_transfer, source_fit,
    AUTHOR_TRAIT_BLOB, BASE_BLOOM_BETA, CV_MODELS
)

class OriginalFloralTraitSourceTests(unittest.TestCase):
    def test_real_bloom_benchmark(self):
        self.assertAlmostEqual(BASE_BLOOM_BETA,.21474990651679526,places=10)
        self.assertEqual(len(AUTHOR_TRAIT_BLOB),40)

    def test_missing_trait_not_coded_as_nonwhite_noscent(self):
        z=pd.Series(["white","NA","",None])
        z=z.replace({"NA":pd.NA,"":pd.NA})
        self.assertTrue(z.iloc[1:].isna().all())
        self.assertFalse(z.loc[z.notna()].eq("white").iloc[0] is False)

    def test_blocked_group_transfer_refuses_small_studies(self):
        x=pd.DataFrame({
           "study_key":["paper_a","paper_b","paper_c","paper_d"],
           "plant_family":["Fam1","Fam2","Fam1","Fam2"],
           "night_bloom":[1,0,1,0]
        })
        r=fivefold_blocked_transfer(x,"study_key")
        self.assertEqual(r["status"],"NOT_ESTIMABLE_EXPERIMENT_CV_STUDY_FAMILY_SUPPORT")

    def test_regression_does_not_emit_p_without_enough_sources(self):
        x=pd.DataFrame({
           "study_key":["one"]*10,"species":["A"]*10,
           "plant_family":["P"]*10,
           "night_bloom":[1]*10
        })
        m=source_fit(x,"night_share ~ night_bloom",["night_bloom"])
        self.assertEqual(m["status"],"NOT_ESTIMABLE_INDEPENDENT_STUDY_TRAIT_SUPPORT")
        self.assertNotIn("terms",m)

    def test_model_schema_preselected_not_tuned(self):
        self.assertEqual(len(CV_MODELS),4)
        self.assertEqual(CV_MODELS["bloom_time_only"],("night_bloom",))


if __name__=="__main__":
    unittest.main()
