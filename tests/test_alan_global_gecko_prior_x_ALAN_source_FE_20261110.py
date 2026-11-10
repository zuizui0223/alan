"""Synthetic-only reference controls for source-observed gecko×light interaction."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_gecko_prior_photo_x_moth_priorPCNL_FE_20261110 import (
    support_group,fit)

class PriorGeckoWorldSourceGates(unittest.TestCase):
    @staticmethod
    def data():
        rng=np.random.default_rng(1207)
        rec=[]
        for island in range(84):
            g=island<42
            for site in range(6):
                light=float(.1+site*.18+rng.normal(0,.025))
                lat=float(8+island*.21+site*.037+rng.normal(0,.02))
                total=int(rng.integers(8,48))
                count=np.log1p(total)
                response=25+2*light+5*g*light+count*1.0+lat*.3+rng.normal(0,.025)
                rec.append({"island":str(island),"year":2021,"latbin":200+island,
                    "lonbin":site,"island_year":f"{island}|2021",
                    "prior_log1p_PCNL":light,
                    "centre_lat":lat,"area_km2":float(30+island*200),
                    "recorded_moth_family_share_pp":response,
                    "log_all_fourfamily_observations":count,
                    "prior_gecko_documented":g,
                    "prior_scink_documented":island%3==0,
                    "prior_Hemidactylus_documented":island%4==0,
                    "prior_Phelsuma_documented":False,
                    "prior_any_reptile_photo_documented":(g or island%3==0)})
        return pd.DataFrame(rec)

    def test_true_island_year_prior_documentation_constancy(self):
        z,s=support_group(self.data())
        self.assertEqual(s["source_prior_gecko_positive_groups"],42)
        self.assertEqual(s["source_gecko_undocumented_groups"],42)
        self.assertEqual(s["original_islands_with_source_insect_sameyear_gradient"],84)

    def test_fitted_moderation_recoverable_in_synthetic_only(self):
        result=fit(self.data(),"prior_gecko_documented")
        self.assertEqual(result["fit"]["status"],
                         "estimated_exploratory_photo_source_composition_not_pollen_service")
        beta=result["fit"]["terms"]["gecko_documentation_x_light"]["beta"]
        self.assertAlmostEqual(beta,5.,delta=.5)
        self.assertTrue(result["source_reference_not_verified_gecko_absent"])

    def test_no_fake_comparator_when_only_reptile_photos_are_gecko(self):
        z=self.data()
        z["prior_any_reptile_photo_documented"]=z["prior_gecko_documented"]
        result=fit(z,"prior_gecko_documented",only_reptile=True)
        self.assertEqual(result["status"],
                         "prior_documentation_interaction_unidentifiable_too_few_islands")

    def test_early_census_year_not_used_for_predator_prior(self):
        z=self.data()
        old=z.copy()
        old["year"]=2013;old["island_year"]=old.island+"|2013"
        mixed=pd.concat([old,z],ignore_index=True)
        valid,s=support_group(mixed)
        self.assertTrue(valid.year.ge(2015).all())

if __name__=="__main__":
    unittest.main()
