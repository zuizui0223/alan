"""Synthetic only: no invented original pixels enter the ecological results."""
import unittest
import pandas as pd
from scripts.analysis import alan_collect_full8942_island_original_roi9_20261108 as m

class Complete8942PhotoSourceAuditTests(unittest.TestCase):
    def test_source_attempts_versus_three_numeric_roi_and_distinct_accepted_taxa(self):
        rows=[]
        for i in range(6):
            rows.append(dict(photo_id=str(i+1),observation_id=f"orig-{i}",
              target_id="2024|a|I",land_id="I",
              wcvp_accepted_plant_name_id=["A","A","B","B","C","C"][i],
              real_numeric_ROI9=(i in [0,2,4])))
        for i in range(3):
            rows.append(dict(photo_id=str(i+7),observation_id=f"orig-{i+6}",
              target_id="2025|b|J",land_id="J",
              wcvp_accepted_plant_name_id="D",
              real_numeric_ROI9=True))
        q=pd.DataFrame(rows)
        x=m.count_source_sites(q)
        self.assertEqual(x["source_photo_cell_years"],2)
        self.assertEqual(x["actual_nine_palette_measured_photos"],6)
        self.assertEqual(x["source_opportunity_site_years_5attempts_3taxa"],1)
        self.assertEqual(x["actual_5attempt_3numericROI9_3observed_accepted_taxa_site_years"],1)
        q.loc[q.observation_id.eq("orig-4"),"real_numeric_ROI9"]=False
        x=m.count_source_sites(q)
        self.assertEqual(x["actual_5attempt_3numericROI9_3observed_accepted_taxa_site_years"],0)

    def test_one_observation_two_images_cannot_be_two_attempts(self):
        q=pd.DataFrame([
          dict(photo_id="1",observation_id="same",target_id="2024|a|I",
               land_id="I",wcvp_accepted_plant_name_id="A",real_numeric_ROI9=True),
          dict(photo_id="2",observation_id="same",target_id="2024|a|I",
               land_id="I",wcvp_accepted_plant_name_id="A",real_numeric_ROI9=True),
        ])
        with self.assertRaisesRegex(ValueError,"pseudo-replicates"):
            m.count_source_sites(q)

    def test_count_unknown_real_photo_fraction_does_not_create_phenotype(self):
        q=pd.DataFrame([
          dict(photo_id="1",observation_id="1",target_id="2024|a|I",
               land_id="I",wcvp_accepted_plant_name_id="A",real_numeric_ROI9=False),
        ])
        x=m.count_source_sites(q)
        self.assertEqual(x["actual_nine_palette_measured_photos"],0)
        self.assertEqual(x["actual_original_ROI9_unknown_photos"],1)

if __name__=="__main__":
    unittest.main()
