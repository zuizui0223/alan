"""Unit tests use artificial records only to check original photo aggregation gates."""
import unittest
import pandas as pd
from scripts.analysis import alan_tier07_collect_full27507_actual_originalROI9 as a

class FullTier07OriginalRGBTests(unittest.TestCase):
    def test_site_strict_needs_five_original_photo_opportunities_and_three_true_species_successes(self):
        df=pd.DataFrame([
         dict(photo_id=str(i),observation_id=f"new{i}",target_id="2024|c|island",
           land_id="island",wcvp_accepted_plant_name_id=["A","B","C","C","D"][i-1],
           measured_roi9=i in (1,2,3)) for i in range(1,6)
        ] + [
         dict(photo_id=str(i),observation_id=f"new{i}",target_id="2025|a|island2",
           land_id="island2",wcvp_accepted_plant_name_id="E",
           measured_roi9=True) for i in range(6,9)
        ])
        x=a.source_and_actual_sites(df)
        self.assertEqual(len(x),2)
        self.assertEqual(x.opportunity_ge5photo_ge3taxa.sum(),1)
        self.assertEqual(x.actual_strict_ge5photo_ge3trueROI_ge3taxa.sum(),1)
        df.loc[df.observation_id.eq("new3"),"measured_roi9"]=False
        y=a.source_and_actual_sites(df)
        self.assertEqual(y.opportunity_ge5photo_ge3taxa.sum(),1)
        self.assertEqual(y.actual_strict_ge5photo_ge3trueROI_ge3taxa.sum(),0)

    def test_same_observation_in_two_photos_never_counts_twice(self):
        data=pd.DataFrame([
          dict(photo_id="900",observation_id="same",target_id="2024|x|L",
               land_id="L",wcvp_accepted_plant_name_id="A",measured_roi9=True),
          dict(photo_id="901",observation_id="same",target_id="2024|x|L",
               land_id="L",wcvp_accepted_plant_name_id="A",measured_roi9=True),
        ])
        with self.assertRaisesRegex(ValueError,"twice"):
            a.source_and_actual_sites(data)

if __name__=="__main__":
    unittest.main()
