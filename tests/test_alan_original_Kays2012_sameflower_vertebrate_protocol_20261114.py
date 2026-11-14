#!/usr/bin/env python3
"""Kays 2012 original published treatment and source-ecological guards."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_Kays2012_sameflower_vertebrate_protocol_20261114 import (
    STUDY_KEY, REFERENCE, SOURCE_KAYS_ARMS, GAMBOA_LAT, GAMBOA_LON,
    haversine_km, two_endpoint_holm, summary_study_sign
)


class OriginalKaysProtocol(unittest.TestCase):
    def test_original_paper_treatment_counts(self):
        self.assertEqual(STUDY_KEY, "Kays_2012")
        self.assertIn("article/download/907/777",REFERENCE)
        self.assertEqual(SOURCE_KAYS_ARMS["day pollination"],(13,1/13))
        self.assertEqual(SOURCE_KAYS_ARMS["night pollination"],(6,.5))
        self.assertEqual(SOURCE_KAYS_ARMS["open pollination"],(6,4/6))
        self.assertEqual(SOURCE_KAYS_ARMS["complete exclusion"],(13,0.))
        I=((4/6)+0-(1/13)-.5)/((4/6)+0+(1/13)+.5)
        self.assertGreater(I,0.)

    def test_primary_gamboa_geometry(self):
        self.assertAlmostEqual(haversine_km(GAMBOA_LAT,GAMBOA_LON,
                                            GAMBOA_LAT,GAMBOA_LON),0)
        self.assertTrue(9 < GAMBOA_LAT < 10)
        self.assertTrue(-80 < GAMBOA_LON < -79)

    def test_study_grain_not_duplicated_by_sites(self):
        z=pd.DataFrame({
            "study_key":["A","A","A","B","B","C"],
            "species":["plant1","plant1","plant2","plant3","plant3","plant4"],
            "two_window_nonadditivity":[-.1,-.2,.2,.1,.3,-.3],
            "open_relative_to_best":[.1,.2,.1,-.1,-.1,.4]
        })
        a=summary_study_sign(z,"two_window_nonadditivity")
        self.assertEqual(a["papers"],3)
        self.assertEqual(a["original_matched_site_year_metric_outcomes"],6)
        self.assertEqual(a["paper_medians_negative"],2)
        self.assertEqual(a["paper_medians_positive"],1)
        b=two_endpoint_holm(z)
        for k in b["nominal"]:
            self.assertGreaterEqual(b["two_new_endpoints_Holm"][k], b["nominal"][k])
        self.assertTrue(b["source_study_omission_sensitivity_not_independent_confirmatory_evidence"])


if __name__=="__main__":
    unittest.main()
