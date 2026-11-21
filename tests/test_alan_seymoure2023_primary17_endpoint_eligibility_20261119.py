import copy
import json
import unittest
from pathlib import Path
from scripts.analysis.alan_seymoure2023_primary17_endpoint_eligibility_20261119 import classify

REG=Path("results/all_species_forest/biotic_interactions/seymoure2023_plant17_primary_endpoint_tiers_20261119.json")
POLLIN={
    "10.1111/gcb.13371",
    "10.1002/ecs2.2550",
    "10.1038/s41598-020-68667-y",
    "10.1016/j.scitotenv.2021.146244",
    "10.1038/s41467-021-24394-0",
    "10.1093/icb/icab010",
}

def mock_source(ledger):
    arr=[]
    for r in ledger["rows"]:
        doi=r["doi"]
        t=(["pollination","herbivory"] if doi=="10.1038/s41598-020-68667-y"
           else ["pollination"] if doi in POLLIN else ["herbivory"])
        arr.append({"doi":doi,"original_2023_review_source_types":t,
                    "original_2023_review_source_interaction_count":len(t)})
    return arr

class Plant17StageEligibilityTests(unittest.TestCase):
    def setUp(self):
        self.reg=json.loads(REG.read_text())
        self.orig=mock_source(self.reg)

    def test_full_primary_tiers_not_a_homogeneous_plant_meta(self):
        d=classify(self.reg,self.orig,strict=True)
        self.assertEqual(d["priority_original_paper_n"],17)
        self.assertEqual(d["priority_original_review_interaction_n"],18)
        self.assertEqual(d["actual_flowering_plant_fruit_seed_papers_n"],3)
        self.assertEqual(d["plant_resources_growth_without_reproductive_fitness_papers_n"],2)
        self.assertEqual(d["other_stage_papers_n"],12)
        self.assertEqual(d["review_guild_vs_real_primary_outcome_matrix"]["pollination"]["flowering_plant_reproduction"],3)
        self.assertEqual(d["review_guild_vs_real_primary_outcome_matrix"]["herbivory"]["flowering_plant_reproduction"],1)
        self.assertEqual(d["effect_size_meta_gate"]["verified_common_same_endpoint_effects_with_original_SE"],0)
        self.assertTrue(d["effect_size_meta_gate"]["no_meta_pooled_mean_or_confidence_interval_fitted"])
        self.assertEqual(sum(d["mutually_exclusive_primary_measured_endpoints_by_paper"].values()),17)

    def test_wrong_plant_endpoint_or_false_precision_rejected(self):
        c=copy.deepcopy(self.reg)
        d=next(x for x in c["rows"] if x["doi"]=="10.1038/s41467-021-24394-0")
        d["measured_living_flowering_plant_reproductive_fitness"]=True
        with self.assertRaisesRegex(ValueError,"mismapped"):
            classify(c,self.orig,strict=True)
        c=copy.deepcopy(self.reg)
        c["rows"][0]["has_verified_independent_light_effect_SE_for_same_comparable_plant_seed_or_fruit_endpoint"]=True
        with self.assertRaisesRegex(ValueError,"source verified"):
            classify(c,self.orig,strict=True)

    def test_missing_and_duplicate_source_paper_rejected(self):
        with self.assertRaisesRegex(ValueError,"not covered"):
            classify(self.reg,self.orig[:-1],strict=True)
        cp=copy.deepcopy(self.reg)
        cp["rows"][0]["doi"]=cp["rows"][1]["doi"]
        with self.assertRaisesRegex(ValueError,"DOI double counted"):
            classify(cp,self.orig,strict=True)

if __name__=="__main__":
    unittest.main()
