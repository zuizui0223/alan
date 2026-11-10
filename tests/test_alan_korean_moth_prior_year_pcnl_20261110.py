"""No invented moth individuals: only synthetic source-year and PCNL-mask unit tests."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

from scripts.analysis.alan_korean_moth_prior_year_pcnl_17year_20261110 import (
    SURVEY_YEARS,SHARDS,years_for_shard,sample_original_raster)

class RealKoreanMothPCNLGates(unittest.TestCase):
    def test_each_actual_survey_year_allocated_once(self):
        sets=[set(years_for_shard(i)) for i in range(SHARDS)]
        self.assertEqual(set.union(*sets),set(SURVEY_YEARS))
        for i in range(SHARDS):
            for j in range(i+1,SHARDS):self.assertFalse(sets[i]&sets[j])
        self.assertEqual(len(SURVEY_YEARS),17)
        self.assertEqual(max(SURVEY_YEARS)-1,2024)
        self.assertEqual(min(SURVEY_YEARS)-1,2006)
    def test_no_missing_PCNL_promoted_to_darkness(self):
        source=np.array([[0,-32768],[5,50]],dtype=np.int16)
        with MemoryFile() as mem:
            with mem.open(driver="GTiff",width=2,height=2,count=1,
                          dtype="int16",crs="EPSG:4326",
                          transform=from_origin(125.,35.,.05,.05),
                          nodata=-32768) as r:
                r.write(source,1)
                sites=pd.DataFrame({"Longitude":[125.025,125.075],
                                    "Latitude":[34.975,34.975]})
                out=sample_original_raster(r,sites)
        self.assertEqual(out[0],0.)
        self.assertTrue(np.isnan(out[1]))
    def test_out_of_range_shard_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"Unknown"):years_for_shard(-1)

if __name__=="__main__":unittest.main()
