"""True original pollen meta strict unit and taxonomic no-moth guard."""
import unittest,random
from scripts.analysis.alan_foldesi2020_original_study_unit_svd import (
 group_summary,within_study_original_plant_pairs,match_lu
)

class TestStudyUnitPollenSource(unittest.TestCase):
    def test_study_equal_aggregation_counts_original_study_not_visitors(self):
        a=[
          {"study":"A","plant":"P","guild":"Diptera",
           "effect_ln_focal_vs_honeybee":0., "pollinator":"f"},
          {"study":"A","plant":"P","guild":"Diptera",
           "effect_ln_focal_vs_honeybee":0., "pollinator":"f2"},
          {"study":"A","plant":"P","guild":"Diptera",
           "effect_ln_focal_vs_honeybee":0., "pollinator":"f3"},
          {"study":"B","plant":"Q","guild":"Diptera",
           "effect_ln_focal_vs_honeybee":1., "pollinator":"g"},
          {"study":"C","plant":"R","guild":"Diptera",
           "effect_ln_focal_vs_honeybee":2., "pollinator":"h"}
        ]
        x=group_summary(a,"Diptera",random.Random(20261118))
        self.assertEqual(x["independent_original_studies"],3)
        self.assertEqual(x["source_original_obs"],5)
        self.assertAlmostEqual(x["source_original_study_means_corrected_lnRatio"],1)

    def test_same_study_fly_and_bee_is_not_fly_and_moth(self):
        a=[
          {"study":"A","plant":"P","guild":"Diptera",
            "effect_ln_focal_vs_honeybee":.5},
          {"study":"A","plant":"P","guild":"Hymenoptera",
            "effect_ln_focal_vs_honeybee":-.2},
          {"study":"A","plant":"Q","guild":"Diptera",
            "effect_ln_focal_vs_honeybee":-.5}
        ]
        x=within_study_original_plant_pairs(a)
        self.assertEqual(x["source_study_plant_guild_pair_count"],1)
        self.assertTrue(x["no_inference_on_moth_pollen_deposition"])
        self.assertTrue(x["study_plant_pairs"][0]["no_valid_SE_due_shared_honeybee_controls"])
        self.assertAlmostEqual(
            x["descriptive_matched_original_study_plant_fly_vs_wild_bee"][
               "matched_original_study_plant_ratio_median"
            ],__import__("math").exp(.7))

    def test_plant_genus_overlap_is_not_functional_effect(self):
        a=[{"plant":"Pedicularis_other"},
           {"plant":"Carum_carvi"},
           {"plant":"Daucus_carota_ssp_sativus"}]
        x=match_lu(a)
        self.assertEqual(x["number_exact_species_overlap"],1)
        self.assertFalse(x["fly_moth_same_plant_SVD_effects_available"])

if __name__=="__main__":
    unittest.main()
