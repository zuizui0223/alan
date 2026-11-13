#!/usr/bin/env python3
"""Unit tests: never infer unmeasured moth pollination from species traits/photos."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_real93m_direct_flower_mothlike_PCNL_20261113 import (
    subset_support,source_site_model,known_site_comparability, LIGHT, KEY
)

class MacroMorphologyTests(unittest.TestCase):
    def test_missing_shape_is_not_negative_in_known_pair(self):
        d=pd.DataFrame({
            "island":["a","a","a","a","b","b"],
            "year":[2024]*6,
            "latbin":[11,11,12,12,1,1],
            "lonbin":[20,20,21,21,40,40],
            "species":["x","z","x","y","x","z"],
            "known_pair":[True,False,True,True,True,False],
            "broad":[True,False,True,False,True,False],
        })
        sites=d[KEY].drop_duplicates().copy()
        sites[LIGHT]=[0,1,0]
        counts=d.loc[d.known_pair].groupby(KEY).size()
        self.assertEqual(int(counts.sum()),4)
        self.assertEqual(d.broad.sum(),3)
        info=subset_support(d,sites,"broad")
        self.assertEqual(info["independent_positive_islands"],2)

    def test_missing_trait_support_fails_closed(self):
        sites=pd.DataFrame({
            "island":["a","a","b","b"],"year":[2024]*4,
            "latbin":[1,1,2,2],"lonbin":[1,2,1,2],
            "all_accepted_recorded_species":[3]*4
        })
        data=sites[KEY].copy()
        data["broad"]=[False,False,True,False]
        x=source_site_model(data,sites,"broad")
        self.assertEqual(x["status"],"not_estimated_insufficient_independent_positive_islands")
        self.assertNotIn("beta_pp",x)

    def test_within_island_year_site_and_light_gate(self):
        sites=pd.DataFrame({
          "island":["a","a","b","b","c","c"],"year":[2024]*6,
          "latbin":[0,0,1,1,2,2],"lonbin":[3,4,3,4,3,4],
          LIGHT:[.1,.8,.2,.21,.0,.1]
        })
        data=pd.concat([sites[KEY],sites[KEY]],ignore_index=True)
        data["known_pair"]=True
        expected=known_site_comparability(data,sites)
        self.assertEqual(set(expected.island),{"a","c"})
        self.assertEqual(len(expected),4)

if __name__=="__main__":
    unittest.main()
