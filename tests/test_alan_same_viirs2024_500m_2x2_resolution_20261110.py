"""Synthetic checks for same-source original VIIRS 2x2 lossless geometry comparison."""
from __future__ import annotations
import unittest
import numpy as np
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from shapely.geometry import box

from scripts.analysis.alan_same_viirs2024_500m_2x2_resolution_20261110 import (
    subset_and_coarsen)

class TrueVIIRSQuartets(unittest.TestCase):
    def runone(self,raw,geom=None):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=100,height=100,count=1,dtype="int16",
                          crs="EPSG:4326",
                          transform=from_origin(-.2083333333333,.2083333333333,1/240,1/240),
                          nodata=-32768) as src:
                src.scales=(10.0,)
                src.offsets=(0.0,)
                src.write(raw.astype(np.int16),1)
                return subset_and_coarsen(src,geom or box(-.25,-.25,.25,.25),
                                          0.,0.,radius_km=5.)
    def test_all_original_dark_has_equal_coarse_fine(self):
        q=self.runone(np.zeros((100,100),dtype=np.int16))
        self.assertEqual([v["threshold_nW_cm2_sr"] for v in q],[0.,1.,3.])
        for v in q:
            self.assertEqual(v["status"],"same_VIIRS_paired_original_pixels")
            self.assertAlmostEqual(v["fine_dark_land_fraction_on_complete_quartets"],1.)
            self.assertAlmostEqual(v["coarse_1km_majority_dark_land_fraction"],1.)
            self.assertAlmostEqual(v["fine_500m_largest_dark_component_fraction_of_dark"],1.)
    def test_checkerboard_fine_and_coarse_are_distinct_without_adding_pixels(self):
        a=np.indices((100,100)).sum(axis=0)%2
        q=self.runone(a*20)
        self.assertGreater(q[0]["fine_dark_land_fraction_on_complete_quartets"],.45)
        self.assertEqual(q[0]["coarse_1km_majority_dark_land_fraction"],0.)
        self.assertEqual(q[0]["coarse_1km_n_majority_dark_blocks"],0)
        self.assertTrue(q[0]["n_complete_land_quartets"]>=16)
    def test_dark_sea_outside_original_land_is_never_included(self):
        raw=np.zeros((100,100),dtype=np.int16)
        raw[:,:50]=10
        q=self.runone(raw,geom=box(-.3,-.3,-.004,.3))
        if q[0]["status"]=="same_VIIRS_paired_original_pixels":
            self.assertAlmostEqual(q[0]["fine_dark_land_fraction_on_complete_quartets"],0.)
            self.assertAlmostEqual(q[0]["coarse_1km_majority_dark_land_fraction"],0.)
    def test_nodata_mask_not_promoted_to_dark(self):
        raw=np.zeros((100,100),dtype=np.int16)
        raw[47:53,47:53]=-32768
        q=self.runone(raw)
        self.assertEqual(q[0]["status"],"same_VIIRS_paired_original_pixels")
        self.assertEqual(q[0]["n_fine_500m_original_land_pixels"],
                         4*q[0]["n_complete_land_quartets"])
        self.assertAlmostEqual(q[0]["fine_dark_land_fraction_on_complete_quartets"],1.)

if __name__=="__main__":unittest.main()
