"""Signed real wetland greenness checks."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(BASE))
spec=importlib.util.spec_from_file_location("signed",BASE/"alan_signed_wetland_ndvi_20261120.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestSignedNDVI(unittest.TestCase):
    def test_three_seasons_early_late_preserved(self):
        rows=[]
        for y in range(2000,2024):
            for day in (111,161,229):
                rows.append(dict(unit_id=33,y=y,day=day,qa=0,
                    ndvi=.8-.01*(y-2000)))
        df=pd.DataFrame(rows)
        a=m.seasonal_annual_ndvi(df)
        v=m.source_eligibility(a,[33])
        self.assertEqual(len(v),1)
        self.assertAlmostEqual(v.iloc[0].ols_ndvi_per_decade,-.1,places=8)
        self.assertLess(v.iloc[0].early_to_late_ndvi,0)
        self.assertAlmostEqual(v.iloc[0].theilsen_ndvi_per_decade,-.1,places=8)

    def test_incomplete_season_rejected(self):
        rows=[]
        for y in range(2000,2024):
            for day in (111,161):
                rows.append(dict(unit_id=33,y=y,day=day,qa=0,ndvi=.5))
        out=m.seasonal_annual_ndvi(pd.DataFrame(rows))
        self.assertTrue(out.empty)

    def test_reproduction_existing_source_17(self):
        for key,path in (("veg",m.VEG),("site",m.SITE),("corine",m.CLC)):
            self.assertEqual(m.git_blob(path),m.PIN[key])
        z=pd.read_csv(m.SITE)
        self.assertEqual(len(z),17)
        d=m.source_eligibility(m.seasonal_annual_ndvi(pd.read_csv(m.VEG)),z.unit_id)
        self.assertEqual(len(d),17)

if __name__=="__main__":
    unittest.main()
