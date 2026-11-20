"""Scientific correctness guards for CORINE status versus true change."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from shapely.geometry import mapping,Polygon,Point
S=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(S))
spec=importlib.util.spec_from_file_location("cha",S/"alan_corine_official_wetland_change_20261120.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def feature(src,dst,box):
    return {"properties":{"CODE_00":src,"CODE_06":dst},
            "geometry":mapping(box)}
class TestCorine(unittest.TestCase):
    def test_pinned_original_17_sites(self):
        p=m.site_source()
        self.assertEqual(len(p),17)
        self.assertTrue(set(p.habitat.astype(str))<=m.WET)

    def test_source_must_provide_true_transition_codes(self):
        p={"properties":{"code_00":"411"}}
        with self.assertRaisesRegex(RuntimeError,"FIELDS"):
            m.extract_transition(p,m.EPOCHS[0])

    def test_polygon_areas_are_clipped_not_whole_polygons(self):
        from shapely.ops import transform
        import pyproj
        circle=Point(0,0).buffer(1000, resolution=24)
        poly=Polygon([(-3000,-3000),(3000,-3000),(3000,3000),(-3000,3000)])
        backwards=pyproj.Transformer.from_crs(3035,4326,always_xy=True).transform
        f=feature("411","211",transform(backwards,poly))
        a=m.summarize_features([f],circle,m.EPOCHS[0])
        self.assertAlmostEqual(a["inland_wetland_loss_ha"],circle.area/10000,places=3)
        self.assertAlmostEqual(a["inland_wetland_gain_ha"],0,places=6)

    def test_no_transition_is_not_a_wetland_loss(self):
        a=m.summarize_features([],Point(0,0).buffer(1000),m.EPOCHS[0])
        self.assertEqual(a["mapped_changes_count"],0)
        self.assertEqual(a["inland_wetland_loss_ha"],0)
        self.assertEqual(a["inland_wetland_gain_ha"],0)

    def test_gain_and_loss_separated(self):
        from shapely.ops import transform
        import pyproj
        circle=Point(0,0).buffer(1000,resolution=24)
        b=pyproj.Transformer.from_crs(3035,4326,always_xy=True).transform
        left=transform(b,Point(-500,0).buffer(100))
        right=transform(b,Point(500,0).buffer(100))
        f=[feature("411","211",left),feature("321","412",right)]
        v=m.summarize_features(f,circle,m.EPOCHS[0])
        self.assertEqual(v["mapped_changes_count"],2)
        self.assertAlmostEqual(v["inland_wetland_loss_ha"],v["inland_wetland_gain_ha"],places=3)
if __name__=="__main__":
    unittest.main()
