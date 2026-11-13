#!/usr/bin/env python3
"""No p-value for a chosen rare island or apparent non-pollinator."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_diel_39island_leaveout_selfing_20261113 import (
    KEY,LIGHT,supported,fit_sample,LEDGER_SHA,BASELINE
)

class DirectFunctionalStabilityTests(unittest.TestCase):
    def setUp(self):
        self.s=pd.DataFrame({
         "island":["x","x","x","z","z","y","y"],
         "year":[2024]*7,
         "latbin":[1,1,1,8,8,4,4],
         "lonbin":[1,2,3,1,2,1,2],
         LIGHT:[0.1,0.3,0.9,1.0,1.2,0.4,0.42],
         "night":[2,1,0,1,3,1,1],
         "day":[1,1,1,2,1,1,0],
         "all_accepted_recorded_species":[3,3,3,4,5,3,3],
         "all_accepted_observation_uuids":[12]*7
        })
    def test_eligibility_is_same_island_year_geographic_light_range(self):
        q,supp=supported(self.s)
        self.assertEqual(supp["positive_sites_both"],5)
        self.assertEqual(supp["eligible_island_years"],2)
        self.assertEqual(supp["evaluated_sites"],4)
        self.assertEqual(supp["evaluated_islands"],2)
        self.assertTrue(q.night.ge(1).all() and q.day.ge(1).all())

    def test_insufficient_independent_islands_is_nonestimable(self):
        x=fit_sample(self.s)
        self.assertEqual(x["status"],"not_estimable_common_source_geography")
        self.assertNotIn("estimate",x)
    def test_predeclared_source_identity(self):
        self.assertEqual(len(LEDGER_SHA),64)
        self.assertAlmostEqual(BASELINE,-1.1747574942083159,places=9)

if __name__=="__main__":
    unittest.main()
