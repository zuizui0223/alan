"""Tests: baseline-water calipers and outcome-blind original matching."""
import importlib.util
import pathlib
import sys
import unittest
import numpy as np
import pandas as pd

DIR=pathlib.Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(DIR))
spec=importlib.util.spec_from_file_location("baseline_match",
    DIR/"alan_strict_baseline_water_20261120.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def synthetic():
    rows=[]
    for block in range(15):
        for j in range(8):
            bright=(j%2==0)
            rows.append(dict(unit_id=block*100+j,block=block,
              light=(2.0+j*.1) if bright else (-2.0-j*.1),
              built=.2+.01*j,built0=2.+.1*j,light0=10.+.1*j,
              start=2000,span=20,temp_trend=.1,rain_trend=.1,
              effort=.02,days=.02,spring=0.,summer=0.,
              freq0=10.+.01*j,area0=10.+.01*j,
              freq_trend=.1,area_trend=.1,
              rate_difference=-2. if bright else 0.,
              candidate_rate=-1. if bright else 0.,
              reference_rate=0.))
    return pd.DataFrame(rows)

class TestStrict(unittest.TestCase):
    def test_no_outcome_enters_matching(self):
        x=synthetic()
        p,_,_,_=m.stricter_pairs(x,.5)
        self.assertGreaterEqual(len(p),30)
        y=x.copy()
        y["rate_difference"]=np.arange(len(y)).astype(float)[::-1]
        y["candidate_rate"]=y.rate_difference
        q,_,_,_=m.stricter_pairs(y,.5)
        pd.testing.assert_frame_equal(
            p[["bright_site_id","comparison_site_id"]].reset_index(drop=True),
            q[["bright_site_id","comparison_site_id"]].reset_index(drop=True))

    def test_strict_baseline_caliper_excludes_different_water(self):
        x=synthetic()
        x.loc[x.light.gt(0),"area0"]+=100
        p,_,_,_=m.stricter_pairs(x,.25)
        self.assertTrue(p.empty)

    def test_caliper_is_fixed_and_water_is_pre_exposure(self):
        x=synthetic()
        p,_,_,limits=m.stricter_pairs(x,.5)
        self.assertLess(limits["area0"],m.maxgap(x.area0,1.))
        self.assertLess(limits["freq0"],m.maxgap(x.freq0,1.))
        with self.assertRaisesRegex(ValueError,"PRESPECIFIED"):
            m.stricter_pairs(x,.75)

    def test_original_187_source_unchanged(self):
        frame,_=m.read_source()
        p,_,_=m.pair_sites(frame,False)
        self.assertEqual(len(p),187)
        self.assertEqual(p.block.nunique(),30)

if __name__=="__main__":
    unittest.main()
