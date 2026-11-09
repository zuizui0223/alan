"""Synthetic contract tests, never substitute for original RGB acquisition."""
import unittest
import pandas as pd
from scripts.acquisition import alan_tier08_24576_balanced_area_light_photo_source_20261109 as m

class Tier08MassivePrecolourTests(unittest.TestCase):
    @staticmethod
    def source():
        rows=[];k=1
        for island,area in (("small","20"),("medium","200"),("large","2000")):
            for site,light in (("x","0.1"),("y","2.0"),("z","1.0")):
                for i in range(8):
                    rows.append({"photo_id":str(k),"observation_id":"uuid-"+str(k),
                        "land_id":island,"target_id":"2024|"+island+site+"|"+island,
                        "year":"2024","wcvp_accepted_plant_name_id":str(i%5),
                        "log_light":light,"area_km2":area})
                    k+=1
        return pd.DataFrame(rows)

    def test_balanced_colour_blind_original_photo_selection_deterministic(self):
        a=self.source()
        s,c=m.frozen_area_light_queue(a,n=48)
        s2,c2=m.frozen_area_light_queue(a.sample(frac=1,random_state=99),n=48)
        self.assertEqual(s.photo_id.tolist(),s2.photo_id.tolist())
        self.assertEqual(c["selected_distinct_original_GSHHG_islands"],3)
        self.assertEqual(c["size"],48)
        self.assertEqual(s.photo_id.nunique(),48)
        self.assertEqual(s.observation_id.nunique(),48)
        self.assertEqual(s.target_id.nunique(),6)
        self.assertEqual(s.land_id.nunique(),3)
        self.assertTrue(s.acquisition_arm.str.startswith("paired_").all())
        self.assertEqual(set(s.area_third),{"small","medium","large"})
        self.assertTrue(s["source_snapshot_date"].eq("2026-09-27").all())

    def test_minimal_sites_and_novel_photo_identifier_coverage(self):
        s,c=m.frozen_area_light_queue(self.source(),n=18)
        self.assertEqual(len(s),18)
        self.assertEqual(int(s.wave.max()),1)
        self.assertEqual(s.photo_id.nunique(),18)
        self.assertTrue(s.shard.eq(0).all())

    def test_rejects_request_exceeding_real_source_opportunities(self):
        with self.assertRaisesRegex(ValueError,"Cannot choose"):
            m.frozen_area_light_queue(self.source(),n=500)

if __name__=="__main__":unittest.main()
