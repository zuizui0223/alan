"""Synthetic contracts for four novel original island photo batches.

No photo images, botanical traits or artificial hues are used in these tests.
"""
from pathlib import Path
import json,tempfile
import unittest

import pandas as pd

from scripts.acquisition import alan_measure_progressive_island_original_roi9_20261108 as a

def fake_schedule():
    rows=[]
    for b in range(1,6):
        for i in range(256):
            pid=str(990000+b*1000+i)
            # Exact signed second batch has 249 distinct GSHHG islands.
            land=(f"{i%249}" if b==2 else f"{i%250}")
            rows.append(dict(
              photo_id=pid,observation_id=f"fake-uuid-{b}-{i}",
              land_id=land,target_id=f"2024|{b}_{i}|{land}",
              source_batch=str(b),colour_measured="False",
              source_arm=a.SOURCE_STUDY_ARM,selection_probability="",
              photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
              photo_license="cc-by",inat_taxon_id=f"{i%20+1}",
              latitude="35.1",longitude="139.2",
              observed_on="2024-07-08",year="2024"))
    return pd.DataFrame(rows)

class OriginalIslandColourBatchTests(unittest.TestCase):
    def test_future_batch_is_256_new_uuid_and_not_first_256(self):
        src=fake_schedule()
        b=a.select_batch(src,2)
        self.assertEqual(len(b),256)
        self.assertEqual(b.land_id.nunique(),249)
        self.assertEqual(b.observation_id.nunique(),256)
        self.assertEqual(b.measurement_id.nunique(),256)
        self.assertEqual(b.shard.value_counts().to_dict(),{0:128,1:128})
        c=a.select_batch(src,3)
        self.assertFalse(set(b.photo_id)&set(c.photo_id))
        self.assertFalse(set(b.observation_id)&set(c.observation_id))
        self.assertFalse(set(b.photo_id)&set(src.loc[src.source_batch.eq("1"),"photo_id"]))
        src.loc[src.source_batch.eq("2"),"photo_id"]=src.loc[src.source_batch.eq("1"),"photo_id"].to_numpy()
        with self.assertRaisesRegex(ValueError,"reused photo"):
            a.select_batch(src,2)

    def test_true_numeric_roi_and_unknown_stay_separate(self):
        import numpy as np
        from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import valid_nine,PALETTE
        frame=pd.DataFrame([
          {**{f"flower_fraction_{p}":(0.65 if p=="white" else 0.35 if p=="pink" else 0) for p in PALETTE},
           "morph":"white","measurement_status":"classified_four_state_morph"},
          {**{f"flower_fraction_{p}":"" for p in PALETTE},
           "morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"},
        ])
        good,fractions=valid_nine(frame)
        self.assertEqual(good.tolist(),[True,False])
        self.assertAlmostEqual(fractions[0].sum(),1.)
        self.assertTrue(np.isnan(fractions[1]).all())

    def test_synthetic_fcp_terminal_receipt_rejoins_all_original_uuids(self):
        q=a.select_batch(fake_schedule(),2)
        with tempfile.TemporaryDirectory() as p:
            root=Path(p)
            source=root/"source"
            (source/"shards").mkdir(parents=True)
            hashes={}
            for part in (0,1):
                s=q.loc[q.shard.eq(part)].copy()
                file=source/"shards"/f"shard-{part:03d}.csv"
                s.to_csv(file,index=False,lineterminator="\n")
                hashes[str(part)]=a.sha(file)
            manifest={"schema":a.SOURCE_SCHEMA,"colour_measured":False,
                      "work_units":2,"worker_ids":[0,1],
                      "photo_shard_sha256":hashes,
                      "prior_original_source_schedule_sha256":a.SCHEDULE_SHA,
                      "original_new_photo_batch_number":2}
            (source/"manifest.json").write_text(json.dumps(manifest))
            measured=root/"measured"
            for part in (0,1):
                s=q.loc[q.shard.eq(part)].copy()
                s["morph"]=["white" if x<80 else "mixed_uncertain" for x in range(128)]
                s["measurement_status"]=s.morph.map({
                     "white":"classified_four_state_morph",
                     "mixed_uncertain":"not_evaluable_roi_or_flip_gate"})
                s["image_sha256"]="f"*64
                s["strict9_white"]=[".8" if x<80 else "" for x in range(128)]
                for name in a.FRACTIONS:
                    s[name]=[".8" if x<80 and name=="flower_fraction_white"
                             else ".2" if x<80 and name=="flower_fraction_pink"
                             else "0" if x<80 else ""
                             for x in range(128)]
                dst=measured/f"shard_{part:03d}"
                dst.mkdir(parents=True)
                s.to_csv(dst/"photos.csv",index=False,lineterminator="\n")
                receipt={"schema":"mois_full_original93m_real_RGB_photo_unit_v1",
                         "shard":part,"attempted":128,
                         "photos_sha256":a.sha(dst/"photos.csv"),
                         "source_sha256":hashes[str(part)],
                         "method":{"measurement_git_blob":a.MODEL_BLOB,
                                   "location_blind":True,
                                   "persisted_image_pixels":False}}
                (dst/"receipt.json").write_text(json.dumps(receipt))
            result=a.collect(source,measured,2,root/"out")
            self.assertEqual(result["original_images_downloaded_and_measured"],256)
            self.assertEqual(result["new_distinct_observation_UUIDs"],256)
            self.assertEqual(result["numeric_original_nine_palette_success"],160)
            self.assertEqual(result["technically_unknown_original_ROI9"],96)
            r=pd.read_csv(root/"out"/"ALAN_original_novel_island_photo_ROI9_batch_02.csv.gz")
            self.assertEqual(len(r),256)
            self.assertEqual(int(r.strict9.sum()),160)
            self.assertTrue(r.loc[~r.strict9,"strict9_white"].isna().all())

if __name__=="__main__":
    unittest.main()
