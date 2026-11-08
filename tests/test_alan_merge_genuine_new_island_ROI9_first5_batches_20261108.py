"""Synthetic source-only statistical counting, not botanical/photo observations."""
import unittest
import pandas as pd
from scripts.analysis import alan_merge_genuine_new_island_ROI9_first5_batches_20261108 as a

class ActualGenuineOriginalColourMergeTests(unittest.TestCase):
    def test_five_attempt_and_three_real_accepted_success_threshold(self):
        q=pd.DataFrame([
          {"target_id":"2019|c|island","observation_id":"1","photo_id":"11",
           "real_numeric_ROI9":True,"wcvp_accepted_plant_name_id":"A","land_id":"island"},
          {"target_id":"2019|c|island","observation_id":"2","photo_id":"12",
           "real_numeric_ROI9":True,"wcvp_accepted_plant_name_id":"B","land_id":"island"},
          {"target_id":"2019|c|island","observation_id":"3","photo_id":"13",
           "real_numeric_ROI9":True,"wcvp_accepted_plant_name_id":"C","land_id":"island"},
          {"target_id":"2019|c|island","observation_id":"4","photo_id":"14",
           "real_numeric_ROI9":False,"wcvp_accepted_plant_name_id":"D","land_id":"island"},
          {"target_id":"2019|c|island","observation_id":"5","photo_id":"15",
           "real_numeric_ROI9":False,"wcvp_accepted_plant_name_id":"E","land_id":"island"},
          {"target_id":"2025|d|island2","observation_id":"6","photo_id":"16",
           "real_numeric_ROI9":True,"wcvp_accepted_plant_name_id":"F","land_id":"island2"},
        ])
        x=a.actual_colour_site_support(q)
        self.assertEqual(x["original_source_photo_cell_years"],2)
        self.assertEqual(x["source_cell_years_with_5plus_attempts_3plus_real_WCVP_taxa"],1)
        self.assertEqual(x["actual_strict_5_attempt_3_success_3_accepted_photo_cell_years"],1)
        self.assertEqual(x["independent_source_islands_represented"],2)
        q.loc[q.photo_id.eq("13"),"real_numeric_ROI9"]=False
        x=a.actual_colour_site_support(q)
        self.assertEqual(x["actual_strict_5_attempt_3_success_3_accepted_photo_cell_years"],0)

    def test_reused_observation_UUIDs_cannot_be_counted_as_source_replication(self):
        d=pd.DataFrame([
          dict(target_id="2024|a|A",observation_id="same",photo_id="91",
               real_numeric_ROI9=True,wcvp_accepted_plant_name_id="F",land_id="A"),
          dict(target_id="2024|a|A",observation_id="same",photo_id="92",
               real_numeric_ROI9=True,wcvp_accepted_plant_name_id="F",land_id="A"),
        ])
        with self.assertRaisesRegex(ValueError,"pseudo-replicated"):
            a.actual_colour_site_support(d)

if __name__=="__main__":
    unittest.main()
