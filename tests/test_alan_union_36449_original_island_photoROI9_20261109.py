"""Synthetic cross-tier constraints, not measured ecological data."""
import unittest
import pandas as pd
from scripts.analysis import alan_union_36449_original_island_photoROI9_20261109 as u

def real_shape_synthetic():
    rows=[]
    for i,(tier,species,success) in enumerate([
        ("tier06","A",True),("tier06","B",False),
        ("tier07","B",True),("tier07","C",True),("tier07","D",False)]):
        rows.append({"photo_id":str(8000+i),"observation_id":f"original-{i}",
          "target_id":"2024|56_84|I","land_id":"I","year":"2024",
          "wcvp_accepted_plant_name_id":species,
          "original_source_tier":tier,"real_numeric_ROI9":success})
    return pd.DataFrame(rows)

class CrossTierMeasuredFlowerColourSiteTests(unittest.TestCase):
    def test_combination_creates_strict_measured_site_neither_tier_has_alone(self):
        d=real_shape_synthetic()
        a=u.site_stats(d.loc[d.original_source_tier.eq("tier06")])
        b=u.site_stats(d.loc[d.original_source_tier.eq("tier07")])
        joint=u.site_stats(d)
        self.assertEqual(int(a.actual_strict_5photo_3success_3taxa.sum()),0)
        self.assertEqual(int(b.actual_strict_5photo_3success_3taxa.sum()),0)
        self.assertEqual(int(joint.actual_strict_5photo_3success_3taxa.sum()),1)
        self.assertEqual(u.summary(d,joint)["source_islands_with_actual_strict_cells"],1)
        self.assertEqual(u.summary(d,joint)["original_9palette_numeric_ROI9"],3)

    def test_genuinely_new_observation_UUID_cannot_repeat_between_source_tiers(self):
        d=real_shape_synthetic()
        d.loc[d.photo_id.eq("8002"),"observation_id"]="original-0"
        with self.assertRaisesRegex(ValueError,"reused"):
            u.site_stats(d)

    def test_same_landmass_with_different_years_never_forms_fake_site(self):
        d=real_shape_synthetic()
        d.loc[d.observation_id.isin(["original-2","original-3","original-4"]),"target_id"]="2025|56_84|I"
        j=u.site_stats(d)
        self.assertEqual(len(j),2)
        self.assertEqual(int(j.actual_strict_5photo_3success_3taxa.sum()),0)

if __name__=="__main__":
    unittest.main()
