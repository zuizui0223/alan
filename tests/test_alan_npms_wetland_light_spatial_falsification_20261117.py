"""Synthetic-only NPMS wetland light-trend and shared site-year falsification checks."""
import unittest

import numpy as np
import pandas as pd
from scripts.analysis import alan_npms_wetland_light_spatial_falsification_20261117 as f


class TestWetlandLightSpatialFalsification(unittest.TestCase):
    def test_source_light_zero_switch_and_trend_sign_without_outcome(self):
        series=[]
        habitat=[]
        for id,light,wet in [
            ("w1",[0.0,1.3,1.4],True),
            ("w2",[1.8,1.6,1.4],True),
            ("o1",[0.9,1.0,1.3],False),
            ("o2",[2.0,1.9,1.8],False),
        ]:
            habitat.append({"location_id":id,"inland_wetland":wet})
            for i,v in enumerate(light):
                series.append({"location_id":id,"year":2017+i,"prior_light":v})
        l=pd.DataFrame(series);h=pd.DataFrame(habitat)
        result=f.plot_light_timeseries(h,l).set_index("location_id")
        self.assertEqual(len(result),4)
        self.assertTrue(bool(result.loc["w1","light_has_zero_and_positive"]))
        self.assertTrue(bool(result.loc["w1","large_annual_light_jump_1logunit"]))
        self.assertFalse(bool(result.loc["w2","light_has_zero_and_positive"]))
        self.assertGreater(result.loc["o1","light_slope_per_year"],0)
        self.assertLess(result.loc["o2","light_slope_per_year"],0)

    def test_shared_source_plotyear_gate_removes_one_guild(self):
        s=pd.DataFrame([
          {"plotyear_method_season":"p1","species":"A","wet7":1.0},
          {"plotyear_method_season":"p1","species":"B","wet7":0.0},
          {"plotyear_method_season":"p2","species":"C","wet7":1.0},
          {"plotyear_method_season":"p2","species":"D","wet7":1.0},
        ])
        z=f.mixed_reselect(s)
        self.assertEqual(len(z),2)
        self.assertEqual(set(z.plotyear_method_season),{"p1"})

    def test_source_support_stop_at_small_wetland_count(self):
        # Does not fabricate estimates from a tiny wetland subset.
        x=pd.DataFrame({
            "location_id":["w","w","o","o"],
            "species":["A","B","C","D"],
            "plotyear_method_season":["wy","wy","oy","oy"],
            "inland_wetland":[True,True,False,False],
            "wet7":[1.,0.,1.,0.],
        })
        out=f.common_support(x)
        self.assertFalse(out["pass"])
        self.assertEqual(out["inland_wetland"]["plots"],1)


if __name__=="__main__":
    unittest.main()
