"""Synthetic safeguards: moth visit evidence is NOT verified effective nocturnal pollination."""
import unittest
import pandas as pd
from scripts.analysis import alan_tier08_moth_positive_source_support_20261109 as m

class MothEvidenceCategories(unittest.TestCase):
    def test_walton_original_study_is_distinguished_from_secondary_citation(self):
        x=pd.DataFrame([
           dict(plant_name="Plantus mothus",guild="moth_family_visit",
                source_family="Noctuidae",relation="visitsFlowersOf",
                study_title="https://github.com/globalbioticinteractions/walton2020/archive/X https://doi.org/10.1098/rsbl.2019.0877"),
           dict(plant_name="Plantus other",guild="moth_family_visit",
                source_family="Sphingidae",relation="pollinates",
                study_title="trophicCH records pollinates"),
           dict(plant_name="Plantus other",guild="butterfly_family_visit",
                source_family="Nymphalidae",relation="visitsFlowersOf",
                study_title="butterfly day source"),
           dict(plant_name="Plantus cited",guild="moth_family_visit",
                source_family="Erebidae",relation="visitsFlowersOf",
                study_title="secondary study said Walton but was not Walton 2020"),
        ])
        e=m.evidence_flags(x).set_index("plant_name")
        self.assertTrue(bool(e.loc["Plantus mothus","walton_direct_nocturnal_pollen_carriage"]))
        self.assertFalse(bool(e.loc["Plantus other","walton_direct_nocturnal_pollen_carriage"]))
        self.assertTrue(bool(e.loc["Plantus other","moth_database_pollinates"]))
        self.assertTrue(bool(e.loc["Plantus other","butterfly_record_positive"]))
        self.assertFalse(bool(e.loc["Plantus cited","walton_direct_nocturnal_pollen_carriage"]))

    def test_true_precolour_site_threshold_does_not_invent_ROI_success(self):
        rows=[]
        for j in range(6):
            rows.append(dict(photo_id=str(j),observation_id="uuid"+str(j),
                 target_id="2024|a|island",land_id="island",
                 wcvp_accepted_plant_name_id=str(j%3),wave=str(1+j//2),
                 log_light=".2",area_third="small",moth_record_positive=True))
        d=pd.DataFrame(rows)
        s=m.support_for(d,"moth_record_positive")
        self.assertEqual(s["precolour_site_years_eligible_for_5_attempt_3_species_minimum"],1)
        self.assertTrue(s["actual_photo_numeric_ROI_success_unknown_until_original_RGB_collector"])
        self.assertEqual(s["source_original_photo_attempts_by_wave"],{"1":2,"2":2,"3":2})

    def test_no_record_never_infers_negative_pollination(self):
        a=pd.DataFrame([
           dict(plant_name="Plantus butterfly",guild="butterfly_family_visit",
                source_family="Pieridae",relation="visitsFlowersOf",
                study_title="study")])
        q=m.evidence_flags(a).iloc[0]
        self.assertFalse(q.moth_record_positive)
        self.assertTrue(q.butterfly_record_positive)
        self.assertFalse(q.walton_direct_nocturnal_pollen_carriage)
if __name__=="__main__":
    unittest.main()
