"""Synthetic regression/geodesic contracts, not a real ecological result."""
import unittest
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import alan_isolation_area_light_20261108 as iso

class IsolationAreaTests(unittest.TestCase):
    def test_mainland_vertex_distance_uses_wrapped_geodesy(self):
        mainland=Point(179.45,0).buffer(.3,resolution=12)
        nearby=Point(-179.55,0).buffer(.20,resolution=12)
        faraway=Point(-165,0).buffer(.2,resolution=12)
        g=gpd.GeoDataFrame({
            "land_id":["100","101","102"],
            "land_class":["mainland","island","island"],
            "area_km2":[9000000,50,60]},
            geometry=[mainland,nearby,faraway],crs=4326)
        q=iso.isolated_mainland_distance_km(g,["101","102"]).set_index("land_id")
        self.assertLess(q.loc["101","distance_mainland_vertex_km"],100)
        self.assertGreater(q.loc["101","distance_mainland_vertex_km"],20)
        self.assertGreater(q.loc["102","distance_mainland_vertex_km"],1000)
        self.assertLess(q.loc["101","distance_mainland_vertex_km"],
                        q.loc["102","distance_mainland_vertex_km"])

    def test_only_outer_coastline_used(self):
        poly=Point(0,0).buffer(1,resolution=12).difference(
            Point(0,0).buffer(.3,resolution=12))
        p=iso.outer_vertices(poly)
        self.assertGreater(len(p),10)
        self.assertTrue(np.all(np.linalg.norm(p,axis=1)>0.9))

    def test_no_missing_distance_or_mainland_misclassification(self):
        d,src=area_join_sample()
        part=d.loc[d.land_id.eq("100")].copy()
        dist=pd.DataFrame({"land_id":["100"],
                           "distance_mainland_vertex_km":[np.nan]})
        with self.assertRaisesRegex(ValueError,"distance join"):
            iso.complete_island_frame(part,dist)
        dist["distance_mainland_vertex_km"]=-1
        with self.assertRaisesRegex(ValueError,"distance join"):
            iso.complete_island_frame(part,dist)

    def test_joint_area_distance_original_within_island_synthetic(self):
        d,area_table=area_join_sample()
        _,within,_=area.join_area_and_select(d,area_table)
        ids=sorted(within.land_id.unique())
        distances=pd.DataFrame({"land_id":ids,
                  "distance_mainland_vertex_km":
                    [1.0+4.0*j for j in [8,3,11,1,9,5,0,10,2,7,4,6]]})
        data=iso.complete_island_frame(within,distances)
        ac=float(np.log10(data.groupby("land_id").area_km2.first()).median())
        dc=float(np.log10(1+data.groupby("land_id")
                           .distance_mainland_vertex_km.first()).median())
        for eq in (False,True):
            a=iso.fit_light_moderation(data,ac,dc,"area_only",eq)
            j=iso.fit_light_moderation(data,ac,dc,"joint",eq)
            self.assertAlmostEqual(
               a["terms"]["log_light:log10_area_c"]["beta"],-1.6,delta=0.8)
            self.assertAlmostEqual(
               j["terms"]["log_light:log10_area_c"]["beta"],-1.6,delta=0.8)
            self.assertTrue(np.isfinite(
               j["terms"]["log_light:log10_distance_c"]["beta"]))
            self.assertEqual(j["n_islands"],12)
            self.assertEqual(j["n_cell_years"],72)
            self.assertEqual(j["weighting"],
                 "equal_total_per_island" if eq else "equal_per_observed_cell_year")

def area_join_sample():
    from importlib.util import spec_from_file_location, module_from_spec
    from pathlib import Path
    p=Path("tests/test_alan_island_area_light_20261108.py")
    spec=spec_from_file_location("original_alan_area_tests",p)
    mod=module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.IslandAreaLightTests.sample()

if __name__=="__main__":
    unittest.main()
