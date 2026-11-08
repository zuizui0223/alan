"""Synthetic verification for ALAN source-only 8942 new-UUID scheduling."""
import unittest
import pandas as pd
from scripts.acquisition import alan_progressive_new_8942_island_original_photo_schedule_20261108 as s

class ProgressiveIslandSourceTests(unittest.TestCase):
    @staticmethod
    def make_source():
        rows=[]
        for i in range(8942):
            pid=str(900000+i)
            land=str(i%369)
            rows.append({
              "photo_id":pid,
              "observation_id":f"fake-uuid-{i:05}",
              "inat_taxon_id":f"taxon-{i%300}",
              "year":"2024",
              "land_id":land,
              "target_id":f"2024|{i}_0|{land}",
              "photo_url_large":
                f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
              "source_arm":s.SOURCE_MASTER_ARM,
              "selection_probability":"",
              "existing_strict_island_cell_years":str((i%369)%5),
            })
        master=pd.DataFrame(rows)
        first=master.iloc[:256].copy()
        first["source_arm"]=s.SOURCE_FIRST_ARM
        return master,first

    def test_schedule_has_no_reused_UUID_and_preserves_first_256(self):
        a,b=self.make_source()
        p=s.schedule(a,b)
        self.assertEqual(len(p),8942)
        self.assertEqual(p.observation_id.nunique(),8942)
        self.assertEqual(p.photo_id.nunique(),8942)
        self.assertEqual(p.land_id.nunique(),369)
        self.assertEqual(p.source_batch.max(),35)
        self.assertEqual(p.groupby("source_batch").size().iloc[-1],238)
        self.assertEqual(p.loc[:255,"observation_id"].tolist(),
                         b.observation_id.tolist())
        self.assertTrue(p.groupby("source_batch").size().iloc[:-1].eq(256).all())
        self.assertEqual(p.loc[p.source_batch.eq(2),"land_id"].nunique(),256)

    def test_all_future_source_batches_are_stable_to_input_order(self):
        a,b=self.make_source()
        x=s.schedule(a,b)
        y=s.schedule(a.sample(frac=1,random_state=9),b)
        self.assertEqual(x.observation_id.tolist(),y.observation_id.tolist())

if __name__=="__main__":
    unittest.main()
