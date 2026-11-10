"""Synthetic-only guards: environmental support never chosen by pollen outcome sign."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    source_stable_group_selection,fit_source)

class OriginalGlobalSourceDarkPatchTests(unittest.TestCase):
    def test_darkness_no_pixel_is_censored_not_zero_contiguity(self):
        src=pd.DataFrame([
           {"island":"a","year":2024,"island_year":"a|2024",
            "threshold_PCNL_native":0.,"status":"observed_satellite_geometry_only",
            "n_dark_raster_pixels":10,"largest_dark_component_fraction_of_dark":.8,
            "dark_land_fraction":.4,"prior_log1p_PCNL":.2,"centre_lat":34.,
            "recorded_moth_family_share_pp":40.,"log_all_fourfamily_observations":2.},
           {"island":"a","year":2024,"island_year":"a|2024",
            "threshold_PCNL_native":0.,"status":"observed_satellite_geometry_only",
            "n_dark_raster_pixels":15,"largest_dark_component_fraction_of_dark":.95,
            "dark_land_fraction":.5,"prior_log1p_PCNL":.5,"centre_lat":34.1,
            "recorded_moth_family_share_pp":42.,"log_all_fourfamily_observations":3.},
           {"island":"a","year":2024,"island_year":"a|2024",
            "threshold_PCNL_native":0.,"status":"observed_satellite_geometry_only",
            "n_dark_raster_pixels":0,"largest_dark_component_fraction_of_dark":np.nan,
            "dark_land_fraction":0.,"prior_log1p_PCNL":.8,"centre_lat":34.2,
            "recorded_moth_family_share_pp":45.,"log_all_fourfamily_observations":4.},
        ])
        selected,support=source_stable_group_selection(src,0.,.05,False)
        self.assertEqual(len(selected),2)
        self.assertEqual(support["retained_islandyear_groups_with_independent_dark_connectivity_contrast"],1)
    def test_not_enough_islands_returns_status_not_coefficient(self):
        x=pd.DataFrame({"island":["one","one"],"island_year":["one|2024"]*2,
              "prior_log1p_PCNL":[0.,.5],"dark_land_fraction":[.8,.3],
              "largest_dark_component_fraction_of_dark":[1.,.4],
              "log_all_fourfamily_observations":[2.,3.],"centre_lat_c":[0.,.3],
              "recorded_moth_family_share_pp":[50.,40.]})
        result=fit_source(x,"recorded_moth_family_share_pp","island_equal",
               ["prior_log1p_PCNL","dark_land_fraction",
                "largest_dark_component_fraction_of_dark",
                "log_all_fourfamily_observations","centre_lat_c"])
        self.assertEqual(result["status"],"insufficient_original_source_islands_for_joint_dark_effect")

if __name__=="__main__":
    unittest.main()
