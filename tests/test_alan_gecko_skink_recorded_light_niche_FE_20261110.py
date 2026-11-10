"""Synthetic-only gecko/skink source composition tests: never an empirical result."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.alan_gecko_skink_recorded_light_niche_FE_20261110 import (
    support_coobserved,true_wls)

class GeckoScinkPhotoLightNicheSafety(unittest.TestCase):
    @staticmethod
    def example():
        rng=np.random.default_rng(86)
        out=[]
        for island in range(55):
            for site in range(6):
                illum=.1+.16*site+rng.normal(0,.012)
                pct=30.+(2.+rng.normal(0,.01))*illum+rng.normal(0,.06)
                n=int(rng.integers(12,30))
                # All numeric fixture values are synthetic test-only data.
                gecko=int(np.clip(round(n*pct/100),1,n-1))
                out.append({"island":str(island),"year":2022,
                            "latbin":100+island,"lonbin":50+site,
                            "island_year":f"{island}|2022",
                            "log1p_prior_PCNL":illum,
                            "gecko_source_photo_count":gecko,
                            "scink_source_photo_count":n-gecko,
                            "original_total_reptile_photo_source_events":n,
                            "gecko_photo_record_percent":100*gecko/n,
                            "log_original_reptile_record_count":np.log1p(n),
                            "centre_lat":20+island*.01+site*.002,
                            "area_km2":100+island,
                           })
        return pd.DataFrame(out)

    def test_genuine_both_photo_categories_required_same_island_year(self):
        x=self.example()
        x.loc[x.island.eq("0"),"gecko_source_photo_count"]=0
        x.loc[x.island.eq("0"),"scink_source_photo_count"]=x.loc[x.island.eq("0"),"original_total_reptile_photo_source_events"]
        z,s=support_coobserved(x)
        self.assertEqual(s["source_retained_independent_islands"],54)
        self.assertFalse(z.island.eq("0").any())
        self.assertTrue(s["no_photo_record_in_cell_means_NOT_DOCUMENTED_not_biological_absence"])

    def test_regression_island_cluster_equalweight(self):
        z,s=support_coobserved(self.example())
        result=true_wls(z)
        self.assertEqual(result["status"],
                         "estimated_exploratory_photo_source_composition_not_pollen_service")
        self.assertEqual(result["source_GSHHG_independent_islands"],55)
        self.assertEqual(result["true_recorded_photo_grid_years"],330)

    def test_no_valid_gradient_cannot_claim_response(self):
        x=self.example()
        x["log1p_prior_PCNL"]=.05
        z,s=support_coobserved(x)
        self.assertTrue(z.empty)
        r=true_wls(z)
        self.assertEqual(r["status"],"not_estimable_empty_authenticated_source")

if __name__=="__main__":
    unittest.main()
