"""Synthetic-only site/year fixed effect projection and source gate tests."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_gecko_scink_repeated_photo_cell_temporalFE_20261110 import (
    repeated_weighted_demean,time_panel,fit_panel)

class OriginalRepeatedReptilePhotoGates(unittest.TestCase):
    def test_two_way_source_FE_projection_removes_photo_location_and_year(self):
        rng=np.random.default_rng(0)
        site=np.repeat(np.arange(55),9)
        year=np.tile(np.arange(2014,2023),55)
        source_obs=rng.normal(size=(len(site),3))
        source_obs+=site[:,None]*.01+year[:,None]*.02
        w=np.ones(len(site))
        centered,iters=repeated_weighted_demean(source_obs,site,year,w)
        self.assertTrue(np.isfinite(centered).all())
        self.assertLess(iters,100)
        self.assertLess(abs(centered.mean(0)).max(),1e-10)
        for j in (0,1,10,54):
            self.assertLess(np.abs(centered[site==j].mean(0)).max(),1e-9)
    def test_original_unrepeated_photo_site_cannot_identify_panel(self):
        frame=pd.DataFrame([
            {"island":str(i),"year":2020,"latbin":30+i,"lonbin":20,
             "log1p_prior_PCNL":i*.01,
             "gecko_source_photo_count":3,
             "scink_source_photo_count":2,
             "gecko_photo_record_percent":60.,
             "original_total_reptile_photo_source_events":5,
             "log_original_reptile_record_count":np.log1p(5),
             "island_year":f"{i}|2020",
             "centre_lat":32.} for i in range(45)])
        z,s=time_panel(frame)
        self.assertTrue(z.empty)
        f=fit_panel(z,s)
        self.assertEqual(f["status"],"unidentifiable_repeated_cell_year_gecko_observation_light")
if __name__=="__main__":
    unittest.main()
