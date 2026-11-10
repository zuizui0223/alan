"""Synthetic-only guards for global same-year moth record / prior dark patch FE.

These tests assert numeric regression and biological gates; NO synthetic
source values enter real analyses or source JSON.
"""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    source_stable_group_selection,fit_source,THRESHOLDS)

class TrueGlobalDarkpatchSourceGates(unittest.TestCase):
    @staticmethod
    def fixture():
        rng=np.random.default_rng(132)
        rows=[]
        for island in range(48):
            for j in range(5):
                light=float(.2+j*.22+rng.uniform(-.02,.02))
                dark=float(.45+.2*rng.uniform(0,1))
                conn=float(.35+.46*rng.uniform(0,1))
                lat=float(35+island*.04+j*.008)
                y=float(18.+2*light+3*dark+5*conn+
                        .5*np.log1p(7+j)+rng.normal(0,.03))
                for threshold in THRESHOLDS:
                    rows.append({
                      "threshold_PCNL_native":threshold,
                      "status":"observed_satellite_geometry_only",
                      "n_dark_raster_pixels":30,
                      "island":str(island),"year":2024,
                      "latbin":int(50+island),"lonbin":int(j),
                      "island_year":f"{island}|2024",
                      "dark_land_fraction":dark,
                      "largest_dark_component_fraction_of_dark":conn,
                      "prior_log1p_PCNL":light,
                      "recorded_moth_family_share_pp":y,
                      "log_all_fourfamily_observations":np.log1p(7+j),
                      "centre_lat":lat,
                    })
        return pd.DataFrame(rows)

    def test_real_comparability_requires_within_island_year_connectivity_range(self):
        d=self.fixture()
        out,diag=source_stable_group_selection(d,0.,.05,False)
        self.assertEqual(diag["retained_actual_independent_islands"],48)
        self.assertEqual(diag["retained_actual_photo_grid_years"],240)
        self.assertTrue((out.n_dark_raster_pixels>0).all())

    def test_same_year_exposure_controls_and_source_photo_bias_separate(self):
        d=self.fixture()
        x,_=source_stable_group_selection(d,0.,.05,False)
        x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
        terms=["prior_log1p_PCNL","dark_land_fraction",
               "largest_dark_component_fraction_of_dark",
               "log_all_fourfamily_observations","centre_lat_c"]
        result=fit_source(x,"recorded_moth_family_share_pp","island_equal",terms)
        self.assertEqual(result["status"],
                         "estimated_exploratory_photo_source_composition_not_pollen_service")
        self.assertEqual(result["source_GSHHG_independent_islands"],48)
        self.assertEqual(result["source_island_year_FEs"],48)
        self.assertEqual(result["true_recorded_photo_grid_years"],240)
        self.assertAlmostEqual(
            result["terms"]["largest_dark_component_fraction_of_dark"]["beta"],5.,delta=.5)

    def test_source_no_detected_dark_does_not_impute_connectedness(self):
        d=self.fixture()
        z=d.index[0]
        d.loc[z,"n_dark_raster_pixels"]=0
        d.loc[z,"largest_dark_component_fraction_of_dark"]=np.nan
        x,_=source_stable_group_selection(d,0.,.05,False)
        self.assertEqual(len(x),239)
        self.assertTrue(x.n_dark_raster_pixels.gt(0).all())

    def test_small_number_islands_does_not_generate_p_values(self):
        d=self.fixture()
        x,_=source_stable_group_selection(d,0.,.05,False)
        x=x.loc[x.island.astype(int)<20].copy()
        x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
        result=fit_source(
            x,"recorded_moth_family_share_pp","island_equal",
            ["prior_log1p_PCNL","dark_land_fraction",
             "largest_dark_component_fraction_of_dark",
             "log_all_fourfamily_observations","centre_lat_c"])
        self.assertEqual(result["status"],
                         "insufficient_original_source_islands_for_joint_dark_effect")

if __name__=="__main__":
    unittest.main()
