#!/usr/bin/env python3
"""30 annual original PCNL sources, exact year chronology and no pseudo estimate."""
import sys,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_authentic_diel_siteyear_PCNL_6shards_20261113 import (
   years,YEARS,SHARDS,START,LAST,N_PAIRS
)
from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
   estimate,LIGHT
)
class OriginalDielTrueYearPCNLTests(unittest.TestCase):
    def test_30_disjoint_original_raster_years(self):
        all_years=[x for i in range(SHARDS) for x in years(i)]
        self.assertEqual(sorted(all_years),list(range(1992,2022)))
        self.assertEqual(len(all_years),len(set(all_years)))
        self.assertTrue(all(len(years(i))==5 for i in range(SHARDS)))
        self.assertEqual((START,LAST,N_PAIRS),(1992,2021,301))
    def test_study_unreplicated_cannot_have_source_p(self):
        d=pd.DataFrame({
          LIGHT:[.2,.3,.4,.5],
          "study_key":["a","b","c","d"],
          "night_share":[.2,.6,.1,.9]
        })
        out=estimate(d,"synthetic","night_share ~ "+LIGHT)
        self.assertEqual(out["status"],"NOT_ESTIMABLE_INDEPENDENT_STUDY_OR_ALAN_RANGE")
        self.assertNotIn("p_nominal_study_cluster",out)
if __name__=="__main__":unittest.main()
