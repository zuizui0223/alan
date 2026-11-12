"""Tests without downloading or pretending to have seen any original image pixels."""
import unittest
import pandas as pd
from scripts.acquisition import alan_tier08_new910_original_rgb_roi9_20261109 as a

class GenuineTier08OriginalPhotoTests(unittest.TestCase):
    def test_sha_and_source_arms_pinned(self):
        self.assertEqual(a.QUEUE_SHA,
                         "9d2e60b9c9bdb2fb56c1c17dd69fca5b6c606d1e64d707020d070dd82e49a552")
        self.assertEqual(a.COUNT,a.REF_COUNT+a.PAIR_COUNT)
        self.assertEqual(a.UNITS,8)
        self.assertEqual(a.MODEL_BLOB,"b1d7091db44682559bbe4b48a9f3b8db648501a6")

    def test_eight_work_shards_and_final_partial_shard(self):
        data=pd.DataFrame({"photo_id":[str(i) for i in range(1,911)],
                           "observation_id":[f"uuid-{i}" for i in range(910)]})
        out=a.assign_shards(data)
        self.assertEqual(out.groupby("shard").size().tolist(),[128]*7+[14])
        self.assertEqual(out.measurement_id.nunique(),910)
        self.assertEqual(out.image_filename.nunique(),910)

    def test_sites_count_only_real_measured_accepted_taxa_not_attempts(self):
        records=[]
        for site in ("A","B"):
            for k in range(5):
                records.append({
                    "target_id":f"2024|{site}|island",
                    "observation_id":f"{site}-{k}",
                    "acquisition_arm":a.PAIR_ARM,
                    "wcvp_accepted_plant_name_id":f"accepted-{k}",
                    "land_id":"island","year":"2024","strict9":k<3,
                })
        t=pd.DataFrame(records)
        result=a.strict_island_sites(t)
        self.assertEqual(result["real_strict_photo_cell_years_targeted_pairs"],2)
        self.assertEqual(result["targeted_island_years_with_two_strict_spatial_sites"],1)

    def test_technical_unknown_is_not_a_measured_species(self):
        records=[]
        for site in ("A","B"):
            for k in range(5):
                records.append({
                    "target_id":f"2024|{site}|island",
                    "observation_id":f"{site}-{k}",
                    "acquisition_arm":a.PAIR_ARM,
                    "wcvp_accepted_plant_name_id":f"accepted-{k}",
                    "land_id":"island","year":"2024",
                    "strict9":k==0,
                })
        result=a.strict_island_sites(pd.DataFrame(records))
        self.assertEqual(result["real_strict_photo_cell_years_targeted_pairs"],0)
        self.assertEqual(result["targeted_island_years_with_two_strict_spatial_sites"],0)

if __name__=="__main__":
    unittest.main()
