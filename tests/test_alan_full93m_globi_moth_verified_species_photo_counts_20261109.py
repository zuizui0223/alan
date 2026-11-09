"""Synthetic unit tests only; source counts in the published results are genuine."""
import unittest
import pandas as pd
from scripts.analysis import alan_full93m_globi_moth_verified_species_photo_counts_20261109 as m

class Source93mGloBIMothCountsTests(unittest.TestCase):
    def test_original_taxon_counts_are_conserved_and_positive_records_not_absence(self):
        count=pd.DataFrame([
            {"inat_taxon_id":"1","species":"A a","n_unique_source_photo_ids":10},
            {"inat_taxon_id":"2","species":"B b","n_unique_source_photo_ids":20},
            {"inat_taxon_id":"3","species":"C c","n_unique_source_photo_ids":30},
            {"inat_taxon_id":"4","species":"D d","n_unique_source_photo_ids":40}])
        wcvp=pd.DataFrame([
            {"inat_taxon_id":"1","wcvp_match_state":"unique_exact_name_to_accepted_id",
             "wcvp_accepted_plant_name_id":"100"},
            {"inat_taxon_id":"2","wcvp_match_state":"unique_exact_name_to_accepted_id",
             "wcvp_accepted_plant_name_id":"200"},
            {"inat_taxon_id":"3","wcvp_match_state":"unique_exact_name_to_accepted_id",
             "wcvp_accepted_plant_name_id":"300"},
            {"inat_taxon_id":"4","wcvp_match_state":"no_accepted_species_match",
             "wcvp_accepted_plant_name_id":""}])
        accepted=pd.DataFrame([
            {"wcvp_accepted_plant_name_id":"100","wcvp_accepted_name":"A a"},
            {"wcvp_accepted_plant_name_id":"200","wcvp_accepted_name":"B b"},
            {"wcvp_accepted_plant_name_id":"300","wcvp_accepted_name":"C c"}])
        globi=pd.DataFrame([
            {"plant_name":"A a","guild":"moth_family_visit","relation":"visitsFlowersOf","study_title":"a moth paper"},
            {"plant_name":"A a","guild":"butterfly_family_visit","relation":"visitsFlowersOf","study_title":"day photo paper"},
            {"plant_name":"B b","guild":"moth_family_visit","relation":"pollinates","study_title":"Walton 2020 https://doi.org/10.1098/rsbl.2019.0877"},
            {"plant_name":"C c","guild":"butterfly_family_visit","relation":"visitsFlowersOf","study_title":"butterfly paper"}])
        report,table=m.compute_counts(count,wcvp,accepted,globi,strict_counts=False)
        self.assertEqual(report["source_photo_count"],100)
        self.assertEqual(report["photo_IDs_in_source_taxa_without_verified_WCVP_accepted_name"],40)
        self.assertEqual(report["evidence_tiers"]["moth_positive"]["original_distinct_photo_IDs"],30)
        self.assertEqual(report["evidence_tiers"]["moth_pollinates_db"]["original_distinct_photo_IDs"],20)
        self.assertEqual(report["evidence_tiers"]["walton2020_night_pollen_carriage"]["original_distinct_photo_IDs"],20)
        self.assertEqual(report["positive_both_moth_and_butterfly_original_photo_IDs"],10)
        self.assertEqual(report["moth_positive_without_recorded_butterfly_original_photo_IDs"],20)
        self.assertEqual(report["evidence_tiers"]["moth_positive"]["512photo_worker_groups_needed_to_measure_all_original_RGB"],1)
        self.assertTrue(report["unknown_GloBI_record_cannot_be_called_moth_absence"])

    def test_no_false_nonflowers_based_on_missing_GloBI_record(self):
        g=pd.DataFrame([{"plant_name":"A a","guild":"butterfly_family_visit",
            "relation":"visitsFlowersOf","study_title":"butterfly"}])
        f=m.evidence(g)
        self.assertEqual(len(f),1)
        self.assertFalse(bool(f.moth_visit_positive.iloc[0]))
        self.assertTrue(bool(f.butterfly_positive.iloc[0]))

    def test_conflicting_WCVP_ID_names_not_silently_multiplied(self):
        count=pd.DataFrame([{"inat_taxon_id":"1","species":"A a","n_unique_source_photo_ids":1}])
        w=pd.DataFrame([{"inat_taxon_id":"1","wcvp_match_state":"unique_exact_name_to_accepted_id",
             "wcvp_accepted_plant_name_id":"1"}])
        a=pd.DataFrame([
            {"wcvp_accepted_plant_name_id":"1","wcvp_accepted_name":"A a"},
            {"wcvp_accepted_plant_name_id":"1","wcvp_accepted_name":"A b"}])
        g=pd.DataFrame([{"plant_name":"A a","guild":"moth_family_visit",
           "relation":"visitsFlowersOf","study_title":"moth"}])
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            m.compute_counts(count,w,a,g,strict_counts=False)

if __name__=="__main__":
    unittest.main()
