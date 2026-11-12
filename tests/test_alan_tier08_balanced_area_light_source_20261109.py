"""Synthetic, pre-colour tests for immutable ALAN tier08 source selection."""
import unittest
import pandas as pd

from scripts.acquisition import alan_tier08_balanced_area_light_source_20261109 as tier

class Tier08OriginalSourceTests(unittest.TestCase):
    def test_source_frame_and_probability_boundary_are_fixed(self):
        self.assertEqual(tier.YEARS, (2019, 2021, 2023, 2024, 2025))
        self.assertEqual(tier.REGISTRY_SHA,
                         "4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b")
        self.assertEqual(tier.REFERENCE_N, 256)
        self.assertEqual(tier.TIER, 8)

    def test_one_photo_per_obs_excludes_prior_and_is_order_independent(self):
        photos = [
            {"photo_id":"10","observation_uuid":"same"},
            {"photo_id":"11","observation_uuid":"same"},
            {"photo_id":"12","observation_uuid":"prior"},
            {"photo_id":"13","observation_uuid":"another"},
        ]
        first = tier.choose_per_observation(photos, {"prior"})
        backwards = tier.choose_per_observation(photos[::-1], {"prior"})
        ids = lambda a: {(x["observation_uuid"],x["photo_id"]) for x in a}
        self.assertEqual(ids(first), ids(backwards))
        self.assertEqual(len(first), 2)
        self.assertNotIn("prior", {r["observation_uuid"] for r in first})

    @staticmethod
    def frame(islands=(("I1", 1.0), ("I2", 100.0), ("I3", 10000.0))):
        records=[]
        for land,area in islands:
            for site,light in (("A",0.1),("B",0.6)):
                for i in range(10):
                    obs=f"{land}_{site}_{i}"
                    records.append({
                        "photo_id":obs, "observation_id":obs,
                        "land_id":land,"year":2024,
                        "target_id":f"2024|{site}|{land}",
                        "wcvp_accepted_plant_name_id":f"taxon{i%5}",
                        "log_light":light,"area_km2":area,
                    })
        return pd.DataFrame(records)

    def test_two_distinct_sameyear_sites_and_balanced_area_bins(self):
        paired,audit=tier.build_paired_queue(self.frame())
        self.assertEqual(audit["candidate_area_light_paired_islands"],3)
        self.assertEqual(audit["selected_paired_islands"],3)
        self.assertEqual(audit["selected_paired_site_years"],6)
        self.assertEqual(len(paired),60)
        self.assertEqual(paired.photo_id.nunique(),60)
        self.assertEqual(paired.observation_id.nunique(),60)
        self.assertTrue(paired.groupby("land_id").site_id.nunique().eq(2).all())
        self.assertEqual(set(paired.area_third),{"small","medium","large"})

    def test_year_repetition_cannot_substitute_for_two_spatial_sites(self):
        m=self.frame((("I1",1.0),))
        # Moving the second site's photographs to 2025 at the same island
        # leaves only one genuine site per island-year.
        mask=m.target_id.eq("2024|B|I1")
        m.loc[mask,"year"]=2025
        m.loc[mask,"target_id"]="2025|B|I1"
        paired,audit=tier.build_paired_queue(m)
        self.assertEqual(len(paired),0)
        self.assertEqual(audit["candidate_area_light_paired_islands"],0)

    def test_absent_within_island_light_gradient_not_fabricated(self):
        m=self.frame((("I1",10.0),))
        m.loc[m.log_light.eq(0.6),"log_light"]=0.11
        paired,audit=tier.build_paired_queue(m)
        self.assertEqual(len(paired),0)
        self.assertEqual(audit["selected_paired_islands"],0)

if __name__=="__main__":
    unittest.main()
