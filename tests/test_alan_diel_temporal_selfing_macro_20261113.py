#!/usr/bin/env python3
"""Prevent temporal pseudo-panel and selfing/moth-study-source double use."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_diel_repeatcell_temporal_20261113 import (
    repeat_panel,longitudinal_regression,LIGHT,KEY
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import (
    portfolio_frame,compatible,weighted_multivariate_year_island_FE
)

class TemporalDirectCompensationSourceTests(unittest.TestCase):
    def test_genuine_revisit_of_same_island_cell_not_two_separate_sites(self):
        x=pd.DataFrame({
          "island":["a","a","a","b","b","b","b"],
          "latbin":[1,1,2,3,3,4,4],
          "lonbin":[10,10,10,10,10,10,10],
          "year":[2023,2024,2024,2023,2024,2023,2024],
          LIGHT:[.20,.24,.4,.10,.25,.1,.4]})
        p,info=repeat_panel(x)
        self.assertEqual(len(p),6)
        self.assertEqual(info["repeat_source_photo_cells"],3)
        self.assertEqual(info["repeated_islands"],2)
        self.assertEqual(info["repeated_cells_with_nonzero_prior_mean_PCNL_change"],3)

    def test_small_panel_no_claim(self):
        x=pd.DataFrame({
          "island":["a","a"],"latbin":[1,1],"lonbin":[1,1],
          "year":[2023,2024],LIGHT:[.10,.22],"log_obs":[2.,2.]})
        s=longitudinal_regression(x,np.array([.4,.5]),"synthetic")
        self.assertEqual(s["status"],"NOT_ESTIMABLE_REPEATED_SOURCE_CELL_SUPPORT")

    def test_experimental_night_species_excluded_from_selfing_denominator(self):
        traits=pd.DataFrame({
            "autonomous_selfing_capacity":["autonomous","absent","autonomous"],
            "self_incompatibility":["SC","SI","SC"]
        },index=["Study sp","Other sp","Third sp"])
        flora=pd.DataFrame({
          "island":["a","a","a","b","b"],
          "latbin":[1,1,1,2,2],"lonbin":[1,1,1,1,1],
          "year":[2024]*5,
          "species":["Study sp","Other sp","Third sp","Study sp","Other sp"],
          "genus":["Study","Other","Third","Study","Other"]})
        sites=flora[KEY].drop_duplicates().copy()
        sites["all_accepted_recorded_species"]=[3,2]
        sites["log_obs"]=[1.,1.]
        sites[LIGHT]=[.3,.4]
        exp=pd.DataFrame({
            "species":["Study sp"],"night_stronger":[True],
            "day_stronger":[False]})
        q,p=portfolio_frame(flora,sites,exp,traits,"autonomous_selfing_capacity")
        self.assertEqual(len(q),2)
        self.assertTrue((q.night_presence==1).all())
        self.assertEqual(q.trait_known.tolist(),[2,1])
        self.assertEqual(q.trait_yes.tolist(),[1,0])
        self.assertEqual(q.trait_no.tolist(),[1,1])
        self.assertNotIn("Study sp",p.species.tolist())

    def test_night_modifier_is_fatal_for_too_few_independent_islands(self):
        q=pd.DataFrame({
          "island":["a","a","b","b"],"year":[2024]*4,
          "night_presence":[0,1,0,1],
          "trait_known":[2]*4,"trait_yes":[1,0,0,1],
          "log_obs":[1]*4, LIGHT:[.2,.4,.3,.5]
        })
        got=weighted_multivariate_year_island_FE(q,q.trait_yes/q.trait_known)
        self.assertEqual(got["status"],"NOT_ESTIMABLE_INDEPENDENT_ISLAND_SOURCE_SUPPORT")

if __name__=="__main__":
    unittest.main()
