"""No invented ecological data: test six-way year partition and raster sampling."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

from scripts.analysis import alan_pcnl_verified_sharded_join_20261109 as s


class TestPCNLSourceJoin(unittest.TestCase):
    def test_33_real_historical_years_partition_without_future_overlap(self):
        yr=sum([s.split_years(i) for i in range(6)],[])
        self.assertEqual(sorted(yr),list(range(1992,2025)))
        self.assertEqual(len(set(yr)),33)
        self.assertEqual(len([y for y in yr if y < 2019]),27)
        self.assertEqual(len([y for y in yr if y < 2025]),33)

    def test_invalid_shard_refuses_source(self):
        with self.assertRaisesRegex(ValueError,"Bad shard"):
            s.split_years(6)

    def test_actual_geotransform_respects_missing_nodata(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"one.tif"
            a=np.array([[0.,2.],[-9999.,4.]],dtype="float32")
            with rasterio.open(f,"w",driver="GTiff",width=2,height=2,count=1,
                dtype="float32",transform=from_origin(0,2,1,1),
                crs="EPSG:4326",nodata=-9999.) as z:z.write(a,1)
            x=pd.DataFrame({"lon":[.5,1.5,.5,1.5],
                            "lat":[1.5,1.5,.5,.5]})
            with rasterio.open(f) as z:
                observed=s.sample_one_raster(z,x)
            self.assertEqual(observed[0],0.)
            self.assertEqual(observed[1],2.)
            self.assertTrue(np.isnan(observed[2]))
            self.assertEqual(observed[3],4.)

    def test_source_provenance_required(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"fake.csv"
            pd.DataFrame({"island":["A"],"year":[2025]}).to_csv(f,index=False)
            with self.assertRaisesRegex(ValueError,"Original full93m source CSV changed"):
                s.source_data(f)

if __name__=="__main__": unittest.main()
