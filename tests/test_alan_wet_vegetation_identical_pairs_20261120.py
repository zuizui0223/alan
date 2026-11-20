"""Tests: original matching must remain identical and NDVI support must be real."""
import importlib.util
import pathlib
import sys
import unittest
import pandas as pd
S=pathlib.Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(S))
spec=importlib.util.spec_from_file_location("wetveg",S/"alan_wet_vegetation_identical_pairs_20261120.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestWetVegetation(unittest.TestCase):
    def test_complete_three_seasons_required(self):
        rows=[]
        for y in range(2000,2024):
            for day in (111,171,231):
                rows.append(dict(unit_id=9,y=y,day=day,qa=0,ndvi=.2+.005*(y-2000)))
        v=m.wetland_vegetation_sites(pd.DataFrame(rows))
        self.assertEqual(len(v),1)
        self.assertEqual(v.iloc[0].valid_years,24)
        self.assertAlmostEqual(v.iloc[0].ndvi_trend,.05,places=7)
        incomplete=pd.DataFrame(rows).query("day == 111")
        self.assertTrue(m.wetland_vegetation_sites(incomplete).empty)

    def test_fixed_pair_missingness_not_imputed_or_rematched(self):
        p=pd.DataFrame([dict(bright_site_id=1,comparison_site_id=2,block=1,outcome_difference=-1),
          dict(bright_site_id=3,comparison_site_id=4,block=2,outcome_difference=-2)])
        v=pd.DataFrame([dict(unit_id=i,ndvi_initial=.4,ndvi_trend=.0) for i in (1,2,3)])
        q,coverage=m.fixed_vegetation_overlap(p,v)
        self.assertEqual(coverage["original_pairs"],2)
        self.assertEqual(coverage["both_arms_ndvi_original_pairs"],1)
        self.assertEqual(m.nested_on_same_pairs(q)["status"],"VEGETATION_SOURCE_SUPPORT_STOP")

    def test_coverage_independent_of_recording_outcomes(self):
        p=pd.DataFrame([dict(bright_site_id=1,comparison_site_id=2,block=1,outcome_difference=-1)])
        v=pd.DataFrame([dict(unit_id=i,ndvi_initial=.4,ndvi_trend=.0) for i in (1,2)])
        _,a=m.fixed_vegetation_overlap(p,v)
        p.outcome_difference=1000
        _,b=m.fixed_vegetation_overlap(p,v)
        self.assertEqual(a,b)

if __name__=="__main__":
    unittest.main()
