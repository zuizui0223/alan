"""Synthetic tests validate GloBI visitor provenance and photographic source contrasts."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_real_oldnew_moth_GloBI_restricted_20261109 as m

class MothPositiveOriginalPhotoStrataTests(unittest.TestCase):
    @staticmethod
    def x():
        old=pd.DataFrame([
          dict(photo_id="1",observation_id="1",wcvp_accepted_plant_name_id="A",
               source_photo_white_fraction=.5,real_numeric_ROI9=True),
          dict(photo_id="2",observation_id="2",wcvp_accepted_plant_name_id="B",
               source_photo_white_fraction=.5,real_numeric_ROI9=True),
          dict(photo_id="3",observation_id="3",wcvp_accepted_plant_name_id="C",
               source_photo_white_fraction=.5,real_numeric_ROI9=True),
          dict(photo_id="4",observation_id="4",wcvp_accepted_plant_name_id=np.nan,
               source_photo_white_fraction=np.nan,real_numeric_ROI9=False)])
        names=pd.DataFrame([
          dict(wcvp_accepted_plant_name_id="A",wcvp_accepted_name="Plantus a"),
          dict(wcvp_accepted_plant_name_id="B",wcvp_accepted_name="Plantus b"),
          dict(wcvp_accepted_plant_name_id="C",wcvp_accepted_name="Plantus c")])
        interactions=pd.DataFrame([
          dict(plant_name="Plantus a",guild="moth_family_visit",
               source_family="Noctuidae",relation="visitsFlowersOf",
               study_title="published moth visit"),
          dict(plant_name="Plantus b",guild="butterfly_family_visit",
               source_family="Nymphalidae",relation="visitsFlowersOf",
               study_title="published butterfly visit"),
          dict(plant_name="Plantus c",guild="moth_family_visit",
               source_family="Sphingidae",relation="pollinates",
               study_title="doi.org/10.1098/rsbl.2019.0877"),
          dict(plant_name="Plantus c",guild="butterfly_family_visit",
               source_family="Nymphalidae",relation="visitsFlowersOf",
               study_title="butterfly publication")])
        return old,names,interactions

    def test_flora_baseline_includes_unknown_taxon_attempts(self):
        x,n,g=self.x()
        q=m.accepted_guilds(x,n,g)
        self.assertEqual(len(q),4)
        self.assertEqual(q.all_original_photographed_flora.sum(),4)
        self.assertEqual(q.moth_record_positive.sum(),2)
        self.assertEqual(q.butterfly_record_positive.sum(),2)
        self.assertEqual(q.moth_record_without_recorded_butterfly.sum(),1)
        self.assertEqual(q.walton_direct_nocturnal_pollen_carriage.sum(),1)
        self.assertTrue(bool(q.loc[q.photo_id.eq("4"),"all_original_photographed_flora"].iloc[0]))
        self.assertFalse(bool(q.loc[q.photo_id.eq("4"),"moth_record_positive"].iloc[0]))

    def test_same_photo_appearing_for_moth_and_butterfly_not_duplicated(self):
        x,n,g=self.x()
        q=m.accepted_guilds(x,n,g)
        self.assertTrue(q.photo_id.is_unique)
        self.assertEqual(int((q.moth_record_positive & q.butterfly_record_positive).sum()),1)

    def test_moth_pollinates_not_same_as_all_moth_positive(self):
        x,n,g=self.x()
        q=m.accepted_guilds(x,n,g)
        self.assertEqual(q.moth_record_positive.sum(),2)
        self.assertEqual(q.moth_database_pollinates.sum(),1)

if __name__=="__main__":unittest.main()
