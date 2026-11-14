#!/usr/bin/env python3
"""Scientific source-value/CI gates; not an invented manuscript raw-data analysis."""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_Grognuz2026_preprint_ALAN_slug_visitor_seed_evidence_gates_20261114 import (
    FILE,SCHEMA,tested_sign,path
)

ROOT=Path(__file__).resolve().parents[1]

class OriginalPreprintPathAudit(unittest.TestCase):
    def get(self):
        return json.loads((ROOT/FILE).read_text())

    def test_both_distinct_damage_pathways_need_full_evidence(self):
        result=path(self.get())
        self.assertEqual(result["status"],
                         "ACTUAL_PREPRINT_REPORTED_RESULT_PATH_AUDIT_NO_RAW_DATA_REANALYSIS")
        self.assertEqual(result["source_version_and_inference_tier"]["independent_publications"],1)
        self.assertFalse(result["source_foliar_damage_mediation_requirements"]
                         ["complete_positive_pollination_mediation_established"])
        self.assertFalse(result["source_floral_damage_mediation_requirements"]
                         ["complete_positive_pollination_mediation_established"])
        self.assertTrue(result["seed_set"]
                        ["NO_actual_seed_assurance_buffering_by_herbivore_state_demonstrated"])
        self.assertTrue(result["source_quality_flags"]
                        ["different_damage_types_must_not_be_merged_in_a_single_mediator"])

    def test_source_sign_check_not_significance_hacking(self):
        self.assertEqual(tested_sign({"beta":.1,"CI95":[-.2,.3]}),
                         "INCONCLUSIVE_REPORTED_CI_CONTAINS_ZERO")
        self.assertEqual(tested_sign({"beta":-1.,"CI95":[-2.,-.3]}),
                         "NEGATIVE_REPORTED_CI_EXCLUDES_ZERO")
        with self.assertRaisesRegex(ValueError,"beta outside"):
            tested_sign({"beta":3.,"CI95":[-1.,1.]})

    def test_cannot_falsely_claim_light_reduced_florivory(self):
        s=self.get()
        s["published_reported_coefficients"]["ALAN_on_florivory_incidence"]
         ["inferential_status"]="LIGHT_REDUCES_FLOWER_DAMAGE"
        with self.assertRaisesRegex(ValueError,"lacks a verified completed"):
            path(s)

if __name__=="__main__":
    unittest.main()
