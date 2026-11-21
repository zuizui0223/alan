import copy
import json
import unittest
from pathlib import Path

from scripts.analysis.alan_sanders23_original_stageA_measurement_process_triage_20261119 import (
    audit, SOURCE
)


class Sanders23StageAProcessIndependenceTests(unittest.TestCase):
    def setUp(self):
        self.reg=json.loads(Path(SOURCE).read_text())
        self.author=[{"source_published_DOIs":r["original_source_DOIs"]}
                     for r in self.reg["source_study_entries"]]

    def test_author_original_23_studies_are_mixed_processes_not_a_meta_sample(self):
        r=audit(self.reg,self.author,strict=True)
        self.assertEqual(r["original_broad_stageA_Keep_StudyID_count"],23)
        self.assertEqual(r["original_study_with_source_DOI_identity_conflict"],1)
        self.assertEqual(r["single_original_stream_drift_preliminary_paper_count"],1)
        self.assertEqual(sum(r["preliminary_outcome_family_counts_not_meta_sample"].values()),23)
        self.assertIsNone(r["total_arthropod_predator_abundance_meta_ready_k"])
        self.assertIsNone(r["riparian_terrestrial_modification_meta_ready_k"])
        self.assertTrue(r["aquatic_mobile_behaviour_is_not_emergent_insect_abundance"])
        self.assertTrue(r["published_author_Source_StudyID_is_not_independent_light_programme"])

    def test_author_bibliography_doi_drift_refused(self):
        author=copy.deepcopy(self.author)
        author[0]["source_published_DOIs"]=["10.0000/fabricated"]
        with self.assertRaisesRegex(ValueError,"DOI list"):
            audit(self.reg,author)

    def test_title_screen_cannot_be_promoted_to_prisma_fulltext(self):
        register=copy.deepcopy(self.reg)
        register["source_study_entries"][0]["final_PRISMA_fulltext_decision"]="INCLUDED"
        with self.assertRaisesRegex(ValueError,"final systematic inclusion"):
            audit(register,self.author)

    def test_source_outcomes_and_experimental_programme_not_auto_retyped(self):
        register=copy.deepcopy(self.reg)
        register["source_study_entries"][3]["preliminary_outcome_family_for_screen_only"]="aquatic_drift_source_NOT_equivalent_to_emergence"
        with self.assertRaisesRegex(ValueError,"classification changed"):
            audit(register,self.author)
        register=copy.deepcopy(self.reg)
        register["source_study_entries"][0]["primary_illumination_programme_id_verified"]=True
        with self.assertRaisesRegex(ValueError,"identity must be"):
            audit(register,self.author)
        register=copy.deepcopy(self.reg)
        register["strict_aquatic_emergence_moderator_meta_eligible_k"]=23
        with self.assertRaisesRegex(ValueError,"remains unknown"):
            audit(register,self.author)


if __name__=="__main__":
    unittest.main()
