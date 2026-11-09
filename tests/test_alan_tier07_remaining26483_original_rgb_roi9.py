"""Contracts only: fake IDs for batching, no fake flower-colour ecological evidence."""
import unittest
import pandas as pd
import numpy as np
from scripts.acquisition import alan_tier07_remaining26483_original_rgb_roi9 as a

def synthetic_original_source():
    rows=[]
    for group in range(27):
        n=883 if group==26 else 1024
        for idx in range(n):
            pid=str(10000000+group*1200+idx)
            rows.append({
             "photo_id":pid,
             "observation_id":f"original_uuid_{group}_{idx}",
             "work_group":str(group),
             "within_work_group":str(idx),
             "land_id":str((idx+group)%563),
             "target_id":f"2025|{group}_{idx}|{(idx+group)%563}",
             "year":"2025","observed_on":"2025-07-21",
             "wcvp_accepted_plant_name_id":str(idx%100+1),
             "inat_taxon_id":str(idx%100+1),
             "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
             "source_priority":"remaining_source_completion",
             "source_arm":a.SOURCE_ARM,
             "source_colour_measured":"False","photo_license":"cc-by",
            })
    return pd.DataFrame(rows)

class Tier07BulkPhotoTests(unittest.TestCase):
    def test_fifty_two_disjoint_photo_tasks_exclude_first_1024(self):
        d=synthetic_original_source()
        self.assertEqual(len(d),27507)
        ids=set()
        obs=set()
        for task in range(a.TASK_COUNT):
            q=a.select_task(d,task)
            expected=371 if task==51 else 512
            self.assertEqual(len(q),expected)
            self.assertEqual(q.shard.nunique(),3 if task==51 else 4)
            self.assertEqual(q.photo_id.nunique(),expected)
            self.assertTrue(q.photo_id.astype(str).str.isnumeric().all())
            self.assertFalse(set(q.photo_id)&ids)
            self.assertFalse(set(q.observation_id)&obs)
            ids.update(q.photo_id)
            obs.update(q.observation_id)
        self.assertEqual(len(ids),26483)
        self.assertEqual(len(obs),26483)
        self.assertFalse(ids&set(d.loc[d.work_group.eq("0"),"photo_id"]))
        self.assertEqual(a.task_source_groups(0),(1,0))
        self.assertEqual(a.task_source_groups(51),(26,1))
        with self.assertRaisesRegex(ValueError,"outside 52"):
            a.task_source_groups(52)

    def test_original_9colour_unknown_kept_unknown(self):
        from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import valid_nine,PALETTE
        f=pd.DataFrame([
           {"morph":"red_pink","measurement_status":"classified_four_state_morph",
            **{f"flower_fraction_{p}":(.55 if p=="pink" else .45 if p=="red" else 0)
               for p in PALETTE}},
           {"morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate",
            **{f"flower_fraction_{p}":"" for p in PALETTE}},
        ])
        actual,nine=valid_nine(f)
        self.assertEqual(actual.tolist(),[True,False])
        self.assertAlmostEqual(nine[0].sum(),1.0)
        self.assertTrue(np.isnan(nine[1]).all())

if __name__=="__main__":
    unittest.main()
