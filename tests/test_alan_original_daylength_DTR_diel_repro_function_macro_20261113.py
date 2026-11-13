#!/usr/bin/env python3
"""Reject fake functional climate inference or unknown study environments."""
import sys
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_daylength_DTR_diel_repro_function_macro_20261113 import (
    normal_study, conditional_models, grouped_transfer, BLOB_DAY, BLOB_DTR,
    ENV_KEYS, MODELS, ENV_COLS
)

class OriginalNightBlooWeatherTests(unittest.TestCase):
    def test_two_author_mojibake_study_names_can_be_reconciled_without_fuzzy_joins(self):
        a=pd.Series(["JÃ¼rgens_2014","Larrea-AlcÃ¡zar_2011","Unrelated_2001"])
        self.assertEqual(normal_study(a).tolist(),
                         ["Jürgens_2014","Larrea-Alcázar_2011","Unrelated_2001"])

    def test_join_original_study_site_year_AND_month(self):
        self.assertTrue(set(["study_key","geo_lat5","geo_lon5","year_exp",
                             "month_start","month_end"]).issubset(ENV_KEYS))
        self.assertEqual(set(ENV_COLS),{"DTR","Daylength"})

    def test_no_study_support_no_weather_interaction_pvalue(self):
        x=pd.DataFrame({
             "study_key":["p","p","q","q"],
             "plant_family":["A","B","A","B"],
             "night_bloom":[0,1,0,1],
             "night_share":[.1,.8,.2,.7]
        })
        z=conditional_models(x)
        self.assertEqual(z["status"],"NOT_ESTIMABLE_INDEPENDENT_REPRODUCTIVE_SOURCE_ENVIRONMENT_OVERLAP")
        v=grouped_transfer(x,"plant_family")
        self.assertEqual(v["status"],"NOT_ESTIMABLE_ORIGINAL_STUDY_FAMILY_ENVIRONMENT_CV_SUPPORT")

    def test_original_git_blob_provenance_and_frozen_modifiers(self):
        self.assertEqual(len(BLOB_DAY),40)
        self.assertEqual(len(BLOB_DTR),40)
        self.assertEqual(len(MODELS),5)
        self.assertIn("night_bloom_x_DTR",MODELS["bloom_x_environment"])

if __name__=="__main__":
    unittest.main()
