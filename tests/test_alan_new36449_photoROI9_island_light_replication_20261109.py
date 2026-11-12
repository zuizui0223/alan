"""Synthetic source contracts only, NOT original flower photo evidence."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as m

class IndependentNightlightReplicationTests(unittest.TestCase):
    def test_accepted_species_equal_photographic_whiteness_not_photo_overweight(self):
        # A: 2 actually measured photos; B: 1 measured; unknown photo
        # is still a source attempt but has no white outcome.
        rows=[
          dict(photo_id="1",observation_id="o1",target_id="2024|10_10|100",
               land_id="100",year=2024,wcvp_accepted_plant_name_id="A",
               real_numeric_ROI9=True,source_photo_white_fraction=.9),
          dict(photo_id="2",observation_id="o2",target_id="2024|10_10|100",
               land_id="100",year=2024,wcvp_accepted_plant_name_id="A",
               real_numeric_ROI9=True,source_photo_white_fraction=.7),
          dict(photo_id="3",observation_id="o3",target_id="2024|10_10|100",
               land_id="100",year=2024,wcvp_accepted_plant_name_id="B",
               real_numeric_ROI9=True,source_photo_white_fraction=.2),
          dict(photo_id="4",observation_id="o4",target_id="2024|10_10|100",
               land_id="100",year=2024,wcvp_accepted_plant_name_id="C",
               real_numeric_ROI9=True,source_photo_white_fraction=.1),
          dict(photo_id="5",observation_id="o5",target_id="2024|10_10|100",
               land_id="100",year=2024,wcvp_accepted_plant_name_id="D",
               real_numeric_ROI9=False,source_photo_white_fraction=float("nan"))
        ]
        all_cells,eligible=m.photo_cell_species_equal(pd.DataFrame(rows))
        self.assertEqual(len(eligible),1)
        self.assertEqual(eligible.photo_attempts.iloc[0],5)
        self.assertEqual(eligible.strict9_success.iloc[0],4)
        self.assertEqual(eligible.accepted_species_n.iloc[0],3)
        # Equal-weight photographed-taxa means, not 4-photo weighted mean.
        self.assertAlmostEqual(eligible.white_pp.iloc[0],100*(.8+.2+.1)/3,places=10)
        self.assertNotAlmostEqual(eligible.white_pp.iloc[0],100*(.9+.7+.2+.1)/4)

    def test_no_new_strict_site_without_three_different_successful_taxa(self):
        rows=[]
        for i in range(6):
            rows.append({
              "photo_id":str(i),"observation_id":str(i),
              "target_id":"2023|25_50|100","land_id":"100","year":2023,
              "wcvp_accepted_plant_name_id":["A","A","B","B","C","D"][i],
              "real_numeric_ROI9":i<4,"source_photo_white_fraction":.4 if i<4 else np.nan
            })
        all_cells,strict=m.photo_cell_species_equal(pd.DataFrame(rows))
        self.assertEqual(len(all_cells),1)
        self.assertEqual(len(strict),0)
        self.assertEqual(int(all_cells.accepted_species_n.iloc[0]),2)

    def test_two_original_image_ids_cannot_count_same_observation_twice(self):
        f=pd.DataFrame([
          dict(photo_id="01",observation_id="SAME",target_id="2024|a|9",
               land_id="9",year=2024,wcvp_accepted_plant_name_id="A",
               real_numeric_ROI9=True,source_photo_white_fraction=.4),
          dict(photo_id="02",observation_id="SAME",target_id="2024|a|9",
               land_id="9",year=2024,wcvp_accepted_plant_name_id="B",
               real_numeric_ROI9=True,source_photo_white_fraction=.8)])
        with self.assertRaisesRegex(ValueError,"duplicated"):
            m.photo_cell_species_equal(f)

    def test_true_missing_original_year_environment_fails_closed(self):
        f=pd.DataFrame([dict(
          target_id="2024|10_10|100",land_id="100",year=2024,
          white_pp=30.,photo_attempts=6,strict9_success=4,
          accepted_species_n=3)])
        env=pd.DataFrame([dict(target_id="2024|10_10|100",land_id="100",
          year=2024,land_class="island",environment_complete=False,
          log_light=.1,log_pop=.2,bio5=20.,log_rain=4.,elevation=30.)])
        land=pd.DataFrame([dict(land_id="100",land_class="island",area_km2=400.)])
        dist=pd.DataFrame([dict(land_id="100",distance_mainland_vertex_km=10.)])
        with self.assertRaisesRegex(ValueError,"lack ORIGINAL-year"):
            m.model_supported(f,env,land,dist)

if __name__=="__main__":
    unittest.main()
