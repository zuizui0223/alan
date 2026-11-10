"""Synthetic VIIRS resolution tests: physical scale, thresholds and missing pixels.

No fake synthetic outcome or source result is substituted for the real 2024 COG.
"""
from __future__ import annotations
import unittest
import numpy as np
from shapely.geometry import box
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from scripts.analysis.alan_dark_refugia_local_pcnl_20261110 import dark_refuge_metrics

class TrueFineVIIRSTests(unittest.TestCase):
    def take(self,raw):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=41,height=41,count=1,dtype="float32",
                          crs="EPSG:4326",transform=from_origin(-.1025,.1025,.005,.005),
                          nodata=-9999.) as dst:
                dst.write(raw.astype(np.float32),1)
                return dark_refuge_metrics(dst,box(-.25,-.25,.25,.25),0.,0.,radius_km=5.,
                                           thresholds=(0.,10.,30.),
                                           min_valid_land_pixels=64,min_coverage=.80)
    def test_realistic_fine_pixel_size_gate_and_3_thresholds(self):
        raw=np.ones((41,41))*30.
        raw[:,20:]=0.
        d=self.take(raw)
        self.assertEqual([x["threshold_PCNL_native"] for x in d],[0.,10.,30.])
        self.assertLess(d[0]["dark_land_fraction"],d[2]["dark_land_fraction"])
        self.assertEqual(d[2]["dark_land_fraction"],1.)
        self.assertEqual(d[0]["status"],"observed_satellite_geometry_only")
    def test_fine_nodata_not_treated_as_dark(self):
        a=np.zeros((41,41))
        a[20,20]=-9999.
        d=self.take(a)
        self.assertEqual(d[0]["n_dark_raster_pixels"],d[0]["n_valid_original_PCNL_pixels"])
        self.assertGreater(d[0]["n_original_island_land_pixels"],d[0]["n_valid_original_PCNL_pixels"])
    def test_all_bright_has_undefined_connected_dark(self):
        d=self.take(np.ones((41,41))*100.)
        self.assertIsNone(d[2]["largest_dark_component_fraction_of_dark"])
        self.assertTrue(d[2]["distance_to_dark_censored_at_radius"])

    def test_original_viirs_int16_nodata_is_not_dark(self):
        raw=np.zeros((41,41),dtype=np.int16)
        raw[20,20]=-32768
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=41,height=41,count=1,dtype="int16",
                          crs="EPSG:4326",transform=from_origin(-.1025,.1025,.005,.005),
                          nodata=-32768) as dst:
                dst.write(raw,1)
                d=dark_refuge_metrics(dst,box(-.25,-.25,.25,.25),0.,0.,
                                      radius_km=5.,thresholds=(0.,10.,30.),
                                      min_valid_land_pixels=64,min_coverage=.80)
        self.assertEqual(d[0]["status"],"observed_satellite_geometry_only")
        self.assertEqual(d[0]["n_dark_raster_pixels"],d[0]["n_valid_original_PCNL_pixels"])
        self.assertGreater(d[0]["n_original_island_land_pixels"],d[0]["n_valid_original_PCNL_pixels"])

if __name__=="__main__":
    unittest.main()
