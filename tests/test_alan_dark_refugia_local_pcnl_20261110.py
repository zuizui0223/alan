"""Synthetic-only geometry checks; these are NOT empirical ALAN results."""
from __future__ import annotations
import unittest

import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from shapely.geometry import box

from scripts.analysis.alan_dark_refugia_local_pcnl_20261110 import dark_refuge_metrics


class TestDarkRefugia(unittest.TestCase):
    def calc(self, data, geom=None, threshold=(0.0,), minimum=8):
        geom = box(-.20, -.20, .20, .20) if geom is None else geom
        with MemoryFile() as mem:
            with mem.open(driver="GTiff", width=31, height=31, count=1,
                          dtype="float32", crs="EPSG:4326",
                          transform=from_origin(-.155,.155,.01,.01),
                          nodata=-9999.) as dst:
                dst.write(data.astype("float32"), 1)
                return dark_refuge_metrics(dst, geom, 0., 0., radius_km=5.,
                                           thresholds=threshold,
                                           min_valid_land_pixels=minimum)

    def test_fully_dark_all_connected(self):
        d=self.calc(np.zeros((31,31)))[0]
        self.assertEqual(d["status"],"observed_satellite_geometry_only")
        self.assertAlmostEqual(d["dark_land_fraction"],1.)
        self.assertAlmostEqual(d["largest_dark_component_fraction_of_dark"],1.)
        self.assertEqual(d["n_four_connected_dark_components"],1)
        self.assertIsNotNone(d["nearest_detected_dark_pixel_km"])

    def test_bright_has_censored_dark_distance(self):
        d=self.calc(np.ones((31,31))*10.)[0]
        self.assertEqual(d["n_dark_raster_pixels"],0)
        self.assertEqual(d["dark_land_fraction"],0.)
        self.assertIsNone(d["largest_dark_component_fraction_of_dark"])
        self.assertTrue(d["distance_to_dark_censored_at_radius"])

    def test_diagonal_patches_not_connected(self):
        arr=np.ones((31,31))*10.
        arr[15,15]=0.
        arr[16,16]=0.
        d=self.calc(arr)[0]
        self.assertEqual(d["n_four_connected_dark_components"],2)
        self.assertAlmostEqual(d["largest_dark_component_fraction_of_dark"],.5,places=3)

    def test_ocean_excluded_even_if_zero_radiance(self):
        arr=np.zeros((31,31))
        arr[:,:15]=10.  # land left half is bright; zero ocean must not count dark
        d=self.calc(arr,geom=box(-.20,-.20,-.005,.20))[0]
        self.assertEqual(d["n_dark_raster_pixels"],0)
        self.assertEqual(d["dark_land_fraction"],0.)

    def test_masked_raster_not_automatically_dark(self):
        arr=np.ones((31,31))*10.
        arr[15,15]=-9999.
        d=self.calc(arr)[0]
        self.assertGreater(d["n_original_island_land_pixels"],
                           d["n_valid_original_PCNL_pixels"])
        self.assertEqual(d["n_dark_raster_pixels"],0)

    def test_dark_fraction_monotonic_across_thresholds(self):
        arr=np.ones((31,31))*2.
        arr[:,16:]=.2
        ds=self.calc(arr,threshold=(0.,1.,3.))
        self.assertEqual(len(ds),3)
        self.assertEqual(ds[0]["dark_land_fraction"],0.)
        self.assertGreater(ds[1]["dark_land_fraction"],0.)
        self.assertAlmostEqual(ds[2]["dark_land_fraction"],1.)

    def test_tiny_island_inadequate_support_not_dark(self):
        arr=np.zeros((31,31))
        d=self.calc(arr,geom=box(-.006,-.006,.006,.006),minimum=16)[0]
        self.assertEqual(d["status"],"insufficient_original_raster_land_support")
        self.assertNotIn("dark_land_fraction",d)

    def test_antimeridian_fails_closed(self):
        arr=np.zeros((31,31))
        with self.assertRaisesRegex(ValueError,"Antimeridian"):
            with MemoryFile() as mem:
                with mem.open(driver="GTiff",width=31,height=31,count=1,
                              dtype="float32",crs="EPSG:4326",
                              transform=from_origin(-.155,.155,.01,.01),
                              nodata=-9999.) as dst:
                    dst.write(arr.astype("float32"),1)
                    dark_refuge_metrics(dst,box(-.3,-.3,.3,.3),179.9,0.)


if __name__=="__main__":
    unittest.main()
