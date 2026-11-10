"""Pure geometric and chronology checks; synthetic points are not empirical GBIF."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

from scripts.analysis.alan_global_GBIF_true4family_priorPCNL_20261110 import (
    OBS_YEARS,SHARDS,years,same_island_photo_grid_centers,genuine_values)

class OriginalGlobalGBIFHistoricalLightSourceGates(unittest.TestCase):
    def test_original_13_observation_years_and_four_disjoint_shards(self):
        self.assertEqual(OBS_YEARS,list(range(2013,2026)))
        bands=[set(years(j)) for j in range(SHARDS)]
        self.assertEqual(set.union(*bands),set(OBS_YEARS))
        self.assertTrue(all(not (bands[j]&bands[k]) for j in range(SHARDS)
                            for k in range(j+1,SHARDS)))
        self.assertEqual(min(OBS_YEARS)-1,2012)
        self.assertEqual(max(OBS_YEARS)-1,2024)

    def test_original_masked_source_noData_does_not_mean_dark(self):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=2,height=2,count=1,
                          dtype="int16",crs="EPSG:4326",
                          transform=from_origin(125.,35.,.05,.05),
                          nodata=-32768) as src:
                src.write(np.array([[0,-32768],[1,7]],dtype=np.int16),1)
                sites=pd.DataFrame({"longitude_center":[125.025,125.075],
                                    "latitude_center":[34.975,34.975]})
                vals=genuine_values(src,sites)
        self.assertEqual(vals[0],0.)
        self.assertTrue(np.isnan(vals[1]))

    def test_site_land_identity_rejects_centre_over_sea(self):
        # Patch function only with independently fabricated geometry in
        # unit test: real source worker must use SHA-verified GSHHG.
        import unittest.mock as mock
        land=gpd.GeoDataFrame(
            {"land_id":["island1"],"area_km2":[100.],
             "land_class":["island"]},
            geometry=[box(124.9,34.9,125.01,35.05)],crs=4326)
        frame=pd.DataFrame({"island":["island1","island1"],
              "longitude_center":[125.0,125.1],
              "latitude_center":[35.0,35.0]})
        with mock.patch(
              "scripts.analysis.alan_global_GBIF_true4family_priorPCNL_20261110.load_signed_land",
              return_value=land):
            z=same_island_photo_grid_centers(frame,"ignored-in-synthetic-test")
        self.assertEqual(z.center_inside_same_original_GSHHG_island.tolist(),
                         [True,False])
if __name__=="__main__":
    unittest.main()
