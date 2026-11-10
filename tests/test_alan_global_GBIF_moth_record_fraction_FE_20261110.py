"""Synthetic statistical validation: fixed effects and source event specificity."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_global_GBIF_moth_record_fraction_within_island_FE_20261110 import (
    original_group_demean,frame,transform_weights)

class GenuineGlobalGBIFSourceFEUnitTests(unittest.TestCase):
    def test_island_year_FE_demeans_both_response_and_light(self):
        rng=np.random.default_rng(12)
        g=np.repeat(["a|2019","b|2019","a|2020","b|2020"],9)
        x=rng.normal(size=(len(g),3))
        w=np.ones(len(g))
        out=original_group_demean(x,g,w)
        for group in set(g):
            self.assertTrue(np.allclose(out[np.array(g)==group].mean(0),0.))
    def test_source_moth_zero_not_used_as_true_biological_absence(self):
        d=pd.DataFrame([
            {"island":"a","year":2020,"island_year":"a|2020",
             "latbin":100,"lonbin":100,"prior_log1p_PCNL":0.1,
             "all_four_family_records":20,"moth_family_records":10,
             "butterfly_family_records":10,"recorded_moth_family_share_pp":50.,
             "centre_lat":35.1,"area_km2":100.},
            {"island":"a","year":2020,"island_year":"a|2020",
             "latbin":100,"lonbin":101,"prior_log1p_PCNL":0.5,
             "all_four_family_records":20,"moth_family_records":0,
             "butterfly_family_records":20,"recorded_moth_family_share_pp":0.,
             "centre_lat":35.2,"area_km2":100.},
            {"island":"a","year":2020,"island_year":"a|2020",
             "latbin":100,"lonbin":102,"prior_log1p_PCNL":0.7,
             "all_four_family_records":10,"moth_family_records":5,
             "butterfly_family_records":5,"recorded_moth_family_share_pp":50.,
             "centre_lat":35.3,"area_km2":100.}
        ])
        x,info=frame(d,"recorded_moth_family_share_pp",False,5)
        self.assertEqual(len(x),2)
        self.assertTrue((x.moth_family_records>0).all())
        z,info2=frame(d,"recorded_moth_family_share_pp",True,5)
        self.assertEqual(len(z),3)
        self.assertFalse(info2["source_conditional_on_both_recorded_families"])
    def test_equal_island_weight_not_equal_observation_weights(self):
        d=pd.DataFrame({"island":["a"]*5+["b"]*2,
                        "island_year":["a|2021"]*5+["b|2021"]*2})
        w=transform_weights(d,"island_equal")
        self.assertAlmostEqual(w[:5].sum(),1.)
        self.assertAlmostEqual(w[5:].sum(),1.)

if __name__=="__main__":
    unittest.main()
