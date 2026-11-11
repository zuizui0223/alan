"""Synthetic-only tests: NPMS real calendar year pre-exposure and nodata fail closed."""
from __future__ import annotations

import unittest
import numpy as np
import pandas as pd
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

from scripts.analysis.alan_npms_11yr_true_plot_priorPCNL_source_20261110 import (
    YEAR_RANGE, SHARDS, shard_years, raster_sample)

class NPMSPastLightSourceTests(unittest.TestCase):
    def test_all_true_NPMS_observation_years_disjoint_and_past(self):
        xs=[set(shard_years(j)) for j in range(SHARDS)]
        self.assertEqual(set.union(*xs),set(YEAR_RANGE))
        self.assertTrue(all(not(xs[j]&xs[k]) for j in range(SHARDS)
                            for k in range(j+1,SHARDS)))
        self.assertEqual([min(YEAR_RANGE)-1,max(YEAR_RANGE)-1],[2014,2024])
    def test_PCNl_zero_is_not_masked_unknown(self):
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=2,height=2,count=1,dtype="int16",
                          crs="EPSG:4326",transform=from_origin(-1.,52.,.05,.05),
                          nodata=-32768) as src:
                src.write(np.array([[0,-32768],[1,7]],dtype=np.int16),1)
                z=pd.DataFrame({"fixed_lon":[-.975,-.925],
                                "fixed_lat":[51.975,51.975]})
                x=raster_sample(src,z)
        self.assertEqual(x[0],0.)
        self.assertTrue(np.isnan(x[1]))
    def test_unknown_shard_fail_closed(self):
        with self.assertRaises(ValueError):shard_years(-1)

if __name__=="__main__":
    unittest.main()
