"""Synthetic regressions verify direction/support, never substitute for actual ecology."""
from __future__ import annotations
import unittest
import pandas as pd
from scripts.analysis.alan_dark_refugia_photosource_paircheck_20261110 import (
    compare_pairs,site_cluster_inference)

class TestPhotoMothMatchedProxy(unittest.TestCase):
    def test_dark_connectivity_contrast_retains_true_source_distance(self):
        d=pd.DataFrame([
            dict(island="a",threshold_PCNL_native=0.,n_dark_raster_pixels=20,
                 largest_dark_component_fraction_of_dark=.9,dark_land_fraction=.5,
                 local_mean_PCNL_native=2.,lon=0.,lat=0.,
                 GloBI_positive_source_percent=30.,log_observation_UUIDs=2.,
                 log_recorded_accepted_species=2.),
            dict(island="a",threshold_PCNL_native=0.,n_dark_raster_pixels=20,
                 largest_dark_component_fraction_of_dark=.5,dark_land_fraction=.55,
                 local_mean_PCNL_native=2.,lon=.20,lat=0.,
                 GloBI_positive_source_percent=15.,log_observation_UUIDs=2.5,
                 log_recorded_accepted_species=2.5),
            dict(island="b",threshold_PCNL_native=0.,n_dark_raster_pixels=20,
                 largest_dark_component_fraction_of_dark=.7,dark_land_fraction=.9,
                 local_mean_PCNL_native=2.,lon=5.,lat=5.,
                 GloBI_positive_source_percent=10.,log_observation_UUIDs=2.,
                 log_recorded_accepted_species=2.)
        ])
        p=compare_pairs(d,0.)
        self.assertEqual(len(p),1)
        self.assertEqual(p.iloc[0].diff_photo_GloBI_positive_pp,15.)
        self.assertGreater(p.iloc[0].distance_km,10.)
        self.assertLess(p.iloc[0].diff_log_photo_UUIDs,0)
    def test_island_not_pair_is_inference_unit(self):
        pairs=pd.DataFrame([{"island":str(i),"diff_photo_GloBI_positive_pp":i-2}
                            for i in range(10) for j in range(3)])
        s=site_cluster_inference(pairs,"diff_photo_GloBI_positive_pp")
        self.assertEqual(s["independent_original_islands"],10)
        self.assertEqual(s["source_matched_photo_site_pairs"],30)
        self.assertAlmostEqual(s["mean_island_equal_higher_connected_minus_lower"],2.5)
if __name__=="__main__":unittest.main()
