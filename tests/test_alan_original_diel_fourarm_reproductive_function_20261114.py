#!/usr/bin/env python3
"""Source-locked four-arm field reproductive ecology safety tests."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_diel_fourarm_reproductive_function_20261114 import (
    SOURCE_N,ARM,canonical_int,canonical_key_from_raw,canonical_key_from_preflight,
    original_mean,descriptive_study_sign,original_observational_pcNL
)

class AuthenticFourArmControls(unittest.TestCase):
    def test_original_study_source_contract(self):
        self.assertEqual(SOURCE_N,(301,118,123))
        self.assertEqual(ARM,("day pollination","night pollination",
                "open pollination","complete exclusion","hand pollination"))

    def test_true_same_study_original_site_and_experimental_year_keys(self):
        a=pd.DataFrame([{
         "study_ID":"P_2015","plant_species":"Some  thing",
         "treatment_effectiveness_metric":"fruit set",
         "coordinates_lat":"-12.456789123","coordinates_lon":"109.23",
         "year_start":"2007","year_end":"NA","month_start":"May",
         "month_end":"June","effectiveness_units":"proportion"}])
        b=pd.Series({
         "study_key":"P_2015","species":"Some thing","metric":"fruit set",
         "lat":-12.456789123,"lon":109.23,"year_exp":2007,
         "year_end":np.nan,"month_start":"May","month_end":"June",
         "units":"proportion"})
        self.assertEqual(canonical_key_from_raw(a),canonical_key_from_preflight(b))
        self.assertIsNone(canonical_int(np.nan))

    def test_author_original_weighted_mean_within_one_arm(self):
        s=pd.DataFrame({"effectiveness_value":["2","10"],
                        "sample_size":["3","1"]})
        self.assertAlmostEqual(original_mean(s),4.0)
        s.iloc[0,s.columns.get_loc("effectiveness_value")]="NA"
        self.assertIsNone(original_mean(s))

    def test_experiment_study_support_fails_without_fake_p_value(self):
        z=pd.DataFrame({"study_key":["one","two"],
                        "species":["A","B"],
                        "open_vs_closed":[-.2,.3]})
        a=descriptive_study_sign(z,"open_vs_closed")
        self.assertEqual(a["status"],"NOT_ESTIMABLE_INDEPENDENT_BIOLOGICAL_STUDY_SUPPORT")
        self.assertNotIn("source_study_sign_test_two_sided_p_nominal",a)

    def test_secondary_light_not_inferable_with_one_paper(self):
        s=pd.DataFrame({
          "study_key":["one"]*5,"species":["A"]*5,
          "two_window_nonadditivity":[-.3,-.2,.0,.2,.1],
          "pcnl_prior_mean_log1p":[0,.1,.4,.5,1]})
        r=original_observational_pcNL(s,"two_window_nonadditivity")
        self.assertEqual(r["status"],"NOT_ESTIMABLE_TRUE_STUDY_LIGHT_SUPPORT")
        self.assertNotIn("p_nominal",r)

if __name__=="__main__":
    unittest.main()
