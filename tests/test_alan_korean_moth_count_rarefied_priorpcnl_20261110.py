"""Synthetic mathematical checks only, NOT fake biological Korean moth results."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_korean_moth_count_rarefied_priorpcnl_20261110 import (
    expected_rarefaction,valid_within_island_year)

class KoreanMothObservationalGates(unittest.TestCase):
    def test_finite_sample_rarefaction_exact(self):
        self.assertAlmostEqual(expected_rarefaction([10,10],20),2.)
        self.assertAlmostEqual(expected_rarefaction([20,0,0],20),1.)
        self.assertAlmostEqual(expected_rarefaction([9,1],5),1.5)
        self.assertAlmostEqual(expected_rarefaction([1]*100,20),20.)
        self.assertTrue(np.isnan(expected_rarefaction([5,7],20)))
    def test_no_fake_zero_or_negative_species_abundance(self):
        with self.assertRaises(ValueError):expected_rarefaction([20,-1],20)
        with self.assertRaises(ValueError):expected_rarefaction([5.2,20],20)
    def test_within_island_and_same_year_gradient_required(self):
        d=pd.DataFrame([
            {"island_year":"a|2021","island":"a","sites":"a1",
             "prior_year_log1p_PCNL":0.,"elevation_km":.1,"y":1},
            {"island_year":"a|2021","island":"a","sites":"a2",
             "prior_year_log1p_PCNL":.3,"elevation_km":.2,"y":2},
            {"island_year":"b|2021","island":"b","sites":"b1",
             "prior_year_log1p_PCNL":.2,"elevation_km":.1,"y":2},
            {"island_year":"b|2021","island":"b","sites":"b2",
             "prior_year_log1p_PCNL":.2,"elevation_km":.1,"y":2},
            {"island_year":"c|2020","island":"c","sites":"c1",
             "prior_year_log1p_PCNL":.1,"elevation_km":.1,"y":2},
        ])
        x,s=valid_within_island_year(d,"y")
        self.assertEqual(len(x),2)
        self.assertEqual(s["within_island_year_contrast_islands"],1)
        self.assertEqual(set(x.island),{"a"})

if __name__=="__main__":unittest.main()
