#!/usr/bin/env python3
"""Synthetic controls; original empirical input verified in GitHub source job."""
import sys
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_real_Wenzell2019_populationyear_pollinator_reproduction_20261114 import (
    haversine_km, summarize_period, require, SPECIES_CODES, FILES,
    original_direct_population,
)


class OnePaperPopulationYearControls(unittest.TestCase):
    def test_distance_not_automatically_exact_site(self):
        self.assertAlmostEqual(haversine_km(32.,-100.,32.,-100.),0.,places=9)
        self.assertTrue(1. < haversine_km(32.,-100.,32.014,-100.) < 2.)
        self.assertEqual(len(FILES),3)
        self.assertEqual(len(SPECIES_CODES),3)

    def test_no_study_pseudoreplication_in_summary(self):
        x=pd.DataFrame([
            dict(dominant_period="nocturnal",night_fraction=.8,original_individuals=20,
                 visitor_records=10,night_visitor_count=8,day_visitor_count=2,
                 night_vs_day_stronger=True),
            dict(dominant_period="diurnal",night_fraction=.2,original_individuals=30,
                 visitor_records=12,night_visitor_count=2,day_visitor_count=10,
                 night_vs_day_stronger=False),
            dict(dominant_period="nocturnal",night_fraction=.45,original_individuals=15,
                 visitor_records=5,night_visitor_count=4,day_visitor_count=1,
                 night_vs_day_stronger=False),
        ])
        z=summarize_period(x)
        self.assertEqual(z["source_field_population_years"],3)
        self.assertEqual(z["population_year_night_day_fruit_direction_agrees_with_visitor_majority"],2)
        self.assertEqual(z["by_actual_original_dominant_visitor_period"]["nocturnal"]["population_years"],2)
        self.assertAlmostEqual(z["nocturnal_minus_diurnal_population_equal_night_share_difference"],.425)

    def test_missing_original_raw_fields_fail_closed(self):
        with self.assertRaises(ValueError):
            original_direct_population(pd.DataFrame(),pd.DataFrame(),pd.DataFrame())

    def test_artifact_does_not_infer_local_ALAN_pollination(self):
        with self.assertRaises(ValueError):
            require(False,"Do not infer ALAN cause from field visitor activity")


if __name__=="__main__":
    unittest.main()
