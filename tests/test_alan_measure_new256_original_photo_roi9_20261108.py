"""Synthetic tests, never interpreted as photographed ecological data."""
from pathlib import Path
import hashlib,json,tempfile
import unittest
import numpy as np,pandas as pd
from scripts.acquisition import alan_measure_new256_original_photo_roi9_20261108 as m

class GenuineNewOriginalROI9Tests(unittest.TestCase):
    def test_only_successful_original_measured_nine_fraction_is_numeric(self):
        d=pd.DataFrame([
         dict(morph="white",measurement_status="classified_four_state_morph",
              **{f"flower_fraction_{x}":(".7" if x=="white" else ".3" if x=="yellow" else "0")
                 for x in m.PALETTE}),
         dict(morph="mixed_uncertain",measurement_status="flower_roi_unclassifiable",
              **{f"flower_fraction_{x}":"0.1" for x in m.PALETTE}),
        ])
        success,p=m.valid_nine(d)
        self.assertEqual(success.tolist(),[True,False])
        self.assertAlmostEqual(p[0].sum(),1.)
        self.assertAlmostEqual(p[0,0],.7)
        self.assertTrue(np.isnan(p[1]).all())

    def test_freeze_requires_256_real_unique_observations_and_signed_sha(self):
        rows=[]
        for i in range(256):
            pid=str(123000+i);year=2024
            rows.append(dict(photo_id=pid,observation_id=f"uuid-{i}",
              inat_taxon_id="1",source_arm="alan93m_new_island_tier06_conditional_priority",
              photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
              observed_on="2024-07-13",photo_license="cc-by",
              image_extension="jpg",measurement_id=f"alan93m_{pid}",
              image_filename=f"alan93m_{pid}.jpg",land_id=f"{900+i}",
              target_id=f"2024|{i}_20|{900+i}",year=year,shard=i//128,
              selection_probability="",source_snapshot_date="2026-09-27"))
        with tempfile.TemporaryDirectory() as t:
            t=Path(t)
            queue=t/"queue.csv.gz"
            pd.DataFrame(rows).to_csv(queue,index=False,compression="gzip")
            rec=t/"receipt.json"
            rec.write_text(json.dumps({"schema":m.TIER_SOURCE_SCHEMA,
               "signed_source_tier06":{"tier06_genuinely_indexed_photo_ids":597070},
               "first_256_blind_new_original_photo_queue":{"original_observations":256},
               "source_photo_per_observation_selected_before_colour":True,
               "source_image_pixels_downloaded":0,
               "underlying_source_SHA256":{"93m_original_tier_registry_SHA256":"z"*64}}))
            original=m.SOURCE_QUEUE_SHA
            m.SOURCE_QUEUE_SHA=m.sha(queue)
            try:
                out=t/"frozen"
                result=m.freeze(queue,rec,out)
                self.assertEqual(result["work_units"],2)
                self.assertEqual(result["worker_ids"],[0,1])
                for unit in (0,1):
                    df,_=m.worker.input_rows(out,unit)
                    self.assertEqual(len(df),128)
                    self.assertEqual(df.observation_id.nunique(),128)
                broken=pd.DataFrame(rows[:-1])
                broken.to_csv(t/"incomplete.csv.gz",compression="gzip",index=False)
                m.SOURCE_QUEUE_SHA=m.sha(t/"incomplete.csv.gz")
                with self.assertRaisesRegex(ValueError,"256"):
                    m.validate_queue(t/"incomplete.csv.gz",rec)
            finally:
                m.SOURCE_QUEUE_SHA=original

if __name__=="__main__":
    unittest.main()
