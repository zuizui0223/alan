"""Synthetic solar clock geometry checks: do not substitute for true GBIF data."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.alan_true_GECKO_GBIF_solar_daynight_source_20261110 import (
    noaa_solar_altitude,original_clock)

class OriginalGeckoSolarClockGates(unittest.TestCase):
    def test_equinox_equator_solar_day_vs_night(self):
        utc=pd.to_datetime(["2024-03-20T12:00:00Z","2024-03-20T00:00:00Z"],utc=True)
        e=noaa_solar_altitude([0.,0.],[0.,0.],utc)
        self.assertGreater(e[0],85.)
        self.assertLess(e[1],-85.)

    def test_real_UTC_timezone_offset_is_normalized(self):
        a=original_clock(pd.Series([2024,2024,2024,2024]),
                         pd.Series([3,3,3,3]),pd.Series([20,20,20,20]),
                         pd.Series(["12:00:00Z","21:00:00+09:00",
                                    "12:00:00","26:00:00Z"]))
        self.assertEqual(a.notna().sum(),2)
        self.assertEqual(a.iloc[0],a.iloc[1])
        self.assertTrue(pd.isna(a.iloc[2]))
        self.assertTrue(pd.isna(a.iloc[3]))

    def test_local_solar_day_vs_night_at_japanese_island(self):
        # 03:00 UTC is local noon around 135E; 12 UTC is local night.
        utc=pd.to_datetime(["2024-06-21T03:00:00Z",
                             "2024-06-21T12:00:00Z"],utc=True)
        e=noaa_solar_altitude([35.,35.],[135.,135.],utc)
        self.assertGreater(e[0],60.)
        self.assertLess(e[1],-6.)

if __name__=="__main__":
    unittest.main()
