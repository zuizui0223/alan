"""No source photo RGB or actual ecology is created by this test."""
import unittest
import pandas as pd
from scripts.acquisition import alan_site_densification_from_8942_20261108 as m

def fake_8942_source():
    rows=[]
    for i in range(8942):
        pid=str(800000+i)
        if i<1280:
            land=str(i%369)
            target=f"2024|old_{i}|{land}"
            original_batch=str(i//256+1)
        elif i<2280:
            site=(i-1280)//20
            land=str(site%20)
            target=f"2024|dense_{site}|{land}"
            original_batch="6"
        else:
            land=str(i%369)
            target=f"2024|unjoined_{i}|{land}"
            original_batch="6"
        rows.append(dict(photo_id=pid,observation_id=f"test-uuid-{i}",
          land_id=land,target_id=target,
          wcvp_accepted_plant_name_id=f"species_{i%25}",
          photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
          source_batch=original_batch,colour_measured="False",
          selection_probability="",existing_strict_island_cell_years="0"))
    return pd.DataFrame(rows)

class IslandPhotoSiteDensificationTests(unittest.TestCase):
    def test_new_source_targets_are_exactly_distinct_and_colour_blind(self):
        d=fake_8942_source()
        q,info=m.select_source(d,limit=512)
        self.assertEqual(len(q),512)
        self.assertEqual(q.photo_id.nunique(),512)
        self.assertEqual(q.observation_id.nunique(),512)
        self.assertTrue(q.source_arm.eq(m.SOURCE_ARM).all())
        self.assertTrue(q.selection_probability.eq("").all())
        self.assertFalse(set(q.observation_id)&set(d.loc[d.source_batch.isin(("1","2","3","4","5")),"observation_id"]))
        self.assertTrue(q.groupby("land_id").size().le(m.PER_ISLAND_MAX).all())
        self.assertTrue(q.groupby("target_id").size().le(m.PER_SITE_MAX).all())
        self.assertEqual(q.site_wave_source_batch.value_counts().to_dict(),{1:256,2:256})
        self.assertEqual(info["original_photos_in_first_five_batches_excluded"],1280)
        self.assertGreater(m.number_of_potential_strict_source_sites(
          pd.concat([d.loc[d.source_batch.isin(("1","2","3","4","5"))],q])), 0)

    def test_same_real_cell_year_and_three_true_accepted_taxa_required(self):
        d=pd.DataFrame([
          {"target_id":"2025|a|1","observation_id":"1","wcvp_accepted_plant_name_id":"A"},
          {"target_id":"2025|a|1","observation_id":"2","wcvp_accepted_plant_name_id":"B"},
          {"target_id":"2025|a|1","observation_id":"3","wcvp_accepted_plant_name_id":"C"},
          {"target_id":"2025|a|1","observation_id":"4","wcvp_accepted_plant_name_id":"D"},
          {"target_id":"2025|a|1","observation_id":"5","wcvp_accepted_plant_name_id":"E"},
          {"target_id":"2025|b|2","observation_id":"6","wcvp_accepted_plant_name_id":"A"},
          {"target_id":"2025|b|2","observation_id":"7","wcvp_accepted_plant_name_id":"A"},
          {"target_id":"2025|b|2","observation_id":"8","wcvp_accepted_plant_name_id":"A"},
          {"target_id":"2025|b|2","observation_id":"9","wcvp_accepted_plant_name_id":"B"},
          {"target_id":"2025|b|2","observation_id":"10","wcvp_accepted_plant_name_id":"B"},
        ])
        self.assertEqual(m.number_of_potential_strict_source_sites(d),1)
        self.assertEqual(m.number_of_potential_strict_source_sites(d.head(4)),0)

if __name__=="__main__":
    unittest.main()
