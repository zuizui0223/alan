"""Scientific source independence guards for 2019 peer reviewed ALAN positive response."""
import unittest
from scripts.analysis.alan_Macgregor2019_original_zenodo_and_study_independence_20261114 import (
    PAPERS,validate
)

class MacgregorOriginalEvidenceGate(unittest.TestCase):
    def test_study_and_outcome_independence(self):
        validate(PAPERS)
        self.assertEqual(len(PAPERS),6)
        self.assertEqual(len(set(s["research_group"] for s in PAPERS)),5)
        self.assertFalse(PAPERS[-1]["reproduction_observation_measured"])
        self.assertTrue(PAPERS[0]["reproduction_observation_measured"])

    def test_2019_positive_not_assumed_pollen_mass(self):
        m=PAPERS[0]
        self.assertEqual(m["direction_full_night_vs_unlit"],"positive")
        self.assertEqual(m["seed_count_full_night_vs_unlit"],"not resolved")
        self.assertFalse(m["direct_moth_visitor_action_measured_in_same_flower"])
        self.assertTrue(m["night_day_access_experiment"])

    def test_preprint_not_conflated_with_peer_reviewed(self):
        p=PAPERS[3]
        self.assertIn("preprint",p["evidence_tier"])
        self.assertIn("compound",p["direction_full_night_vs_unlit"])

    def test_tibetan_two_papers_not_independent(self):
        p,j=PAPERS[4:]
        self.assertEqual(p["research_group"],j["research_group"])
        self.assertTrue(p["reproduction_observation_measured"])
        self.assertFalse(j["reproduction_observation_measured"])
        self.assertNotEqual(p["doi"],j["doi"])
        self.assertFalse(p["evidence_tier"].endswith("original_numeric_refit"))

    def test_no_corrupted_source_identifiers(self):
        self.assertEqual(len({x["doi"] for x in PAPERS}),len(PAPERS))
        mutated=[dict(x) for x in PAPERS]
        mutated[-1]["doi"]=mutated[-2]["doi"]
        with self.assertRaises(AssertionError):
            validate(mutated)

if __name__=="__main__":
    unittest.main()
