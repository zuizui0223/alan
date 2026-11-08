"""Synthetic contracts for large original-photo ALAN acquisition; no fabricated ecolabels."""
import unittest
import pandas as pd
import numpy as np
from scripts.acquisition import alan_bulk_7662_original_roi9_source_20261108 as m

class BulkOriginalROI9SourceTests(unittest.TestCase):
    @staticmethod
    def make_schedule():
        rows=[]
        for b in range(1,36):
            n=238 if b==35 else 256
            for i in range(n):
                pid=str(b*100000+i)
                obs=f"novel_{b}_{i}"
                rows.append({
                    "observation_id":obs,"photo_id":pid,
                    "inat_taxon_id":"123",
                    "wcvp_accepted_plant_name_id":"wcvp123",
                    "land_id":f"land{(b+i)%369}",
                    "target_id":f"2024|{b}_{i}|land{(b+i)%369}",
                    "year":"2024","observed_on":"2024-07-08",
                    "photo_license":"cc-by","source_batch":str(b),
                    "photo_order_in_batch":str(i),
                    "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                    "source_arm":m.TAXON_SOURCE_ARM,"colour_measured":"False",
                    "selection_probability":""})
        return pd.DataFrame(rows)

    def test_true_7662_unique_original_UUIDs_and_15_warm_inferences(self):
        d=self.make_schedule()
        self.assertEqual(len(d),8942)
        observed=set(d.loc[d.source_batch.astype(int).le(5),"observation_id"])
        source_ids=[]
        seen_photos=set()
        for g in range(m.GROUPS):
            x=m.select_group(d,g)
            self.assertEqual(len(x),m.group_size(g))
            self.assertFalse(set(x.observation_id)&observed)
            self.assertFalse(set(x.photo_id)&seen_photos)
            seen_photos.update(x.photo_id)
            source_ids.extend(x.observation_id)
            self.assertTrue(x.groupby("shard").size().iloc[:3].eq(128).all())
            if g==14:
                self.assertEqual(len(x),494)
                self.assertEqual(x.groupby("shard").size().iloc[-1],110)
        self.assertEqual(len(source_ids),7662)
        self.assertEqual(len(set(source_ids)),7662)

    def test_first_five_original_batches_never_selected_again(self):
        d=self.make_schedule()
        for group in range(m.GROUPS):
            q=m.select_group(d,group)
            self.assertTrue(q.source_batch.astype(int).ge(6).all())
        self.assertEqual(m.group_batches(0),(6,7))
        self.assertEqual(m.group_batches(14),(34,35))
        with self.assertRaisesRegex(ValueError,"outside 15"):
            m.group_batches(15)

    def test_original_ROI9_unknown_is_never_inferred(self):
        from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import valid_nine,PALETTE
        src=pd.DataFrame([
            {**{f"flower_fraction_{k}":(.7 if k=="white" else .3 if k=="yellow" else 0.)
                for k in PALETTE},
             "morph":"white","measurement_status":"classified_four_state_morph"},
            {**{f"flower_fraction_{k}":"" for k in PALETTE},
             "morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"}])
        good,fractions=valid_nine(src)
        self.assertEqual(good.tolist(),[True,False])
        self.assertAlmostEqual(fractions[0,0],.7)
        self.assertTrue(np.isnan(fractions[1]).all())

if __name__=="__main__":
    unittest.main()
