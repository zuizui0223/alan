import copy
import json
import unittest
from pathlib import Path
from scripts.analysis.alan_predator_role_study_heterogeneity_20261118 import synthesize


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results/all_species_forest/biotic_interactions/predator_role_original_published_programme_ledger_20261118.json"


class PredatorRoleStudyHeterogeneityTests(unittest.TestCase):
    def setUp(self):
        self.registry = json.loads(LEDGER.read_text())

    def test_full_programmes_and_opposing_signs_without_meta(self):
        r = synthesize(self.registry)
        self.assertEqual(r["evidence"]["independent_programmes_unique"], 6)
        self.assertEqual(r["evidence"]["distinct_endpoint_entries"], 9)
        self.assertEqual(r["evidence"]["predator_abundance_programmes"], 5)
        self.assertEqual(r["evidence"]["predator_abundance_published_directions"],
                         {"increase": 3, "not_resolved": 1, "decrease": 1})
        self.assertTrue(r["evidence"]["forest_predator_direction_crosses_increase_null_decrease"])
        self.assertTrue(r["evidence"]["no_common_endpoint_random_effects_meta_fitted"])
        self.assertFalse(r["scope"]["systematic_review"])

    def test_brown_is_not_treatment_contrast_or_independent_nights(self):
        r = synthesize(self.registry)["brown2023_original_published_period_arithmetic"]
        self.assertAlmostEqual(r["during_to_pre_descriptive_ratio"], 11.34/4.96)
        self.assertAlmostEqual(r["post_to_pre_descriptive_ratio"], 3.10/4.96)
        self.assertTrue(r["not_a_control_adjusted_BACI_effect"])
        # Predator and herbivore data are shared original pitfall traps, not two experiments.
        q = synthesize(self.registry)["brown2023_predator_herbivore_source_trophic_capture_reweighting"]
        self.assertAlmostEqual(q["predator_to_herbivore_capture_ratio"]["pre"], 4.96/3.28)
        self.assertAlmostEqual(q["predator_to_herbivore_capture_ratio"]["during"], 11.34/1.45)
        self.assertAlmostEqual(q["predator_to_herbivore_capture_ratio"]["post"], 3.10/2.65)
        self.assertAlmostEqual(q["light_on_to_pre_predator_herbivore_ratio_of_ratios"],
                               (11.34/1.45)/(4.96/3.28))
        self.assertTrue(q["not_a_direct_predation_or_plant_damage_effect"])
        self.assertTrue(q["same_original_trap_level_process_not_two_experiments"])

        self.assertTrue(r["repeated_15_nights_are_not_independent_light_replicates"])

    def test_boom_factorization_and_nonsignificant_net_seeds(self):
        r = synthesize(self.registry)["boom2020_spectral_nursery_algebra"]
        expected = (0.75 * 0.47 * 0.43)/(0.73 * 0.72 * 0.56)
        self.assertAlmostEqual(r["green_to_dark_product_ratio"], expected)
        self.assertAlmostEqual(r["green_to_dark_fruit_factor"] * r["green_to_dark_noninfested_factor"]
                               * r["green_to_dark_conditional_fertilization_factor"], expected)
        self.assertFalse(r["overall_net_seed_decline_statistically_resolved"])
        self.assertEqual(r["overall_seed_output_reported_omnibus_p"], 0.63)

    def test_duplicate_endpoint_and_missing_programme_rejected(self):
        tmp = copy.deepcopy(self.registry)
        tmp["rows"].append(copy.deepcopy(tmp["rows"][0]))
        with self.assertRaisesRegex(ValueError, "Duplicated endpoint"):
            synthesize(tmp)
        tmp = copy.deepcopy(self.registry)
        for i, row in enumerate(tmp["rows"]):
            row["programme_id"] = f"fake_independent_{i}"
        with self.assertRaisesRegex(ValueError, "Missing same-programme"):
            synthesize(tmp)

    def test_do_not_invent_exact_newzealand_or_green_significance(self):
        tmp = copy.deepcopy(self.registry)
        row = next(x for x in tmp["rows"] if x["programme_id"] == "Farnworth2026")
        row["effect"]["value"] = .333333333333
        with self.assertRaisesRegex(ValueError, "approximate"):
            synthesize(tmp)
        tmp = copy.deepcopy(self.registry)
        row = next(x for x in tmp["rows"] if x["programme_id"] == "Boom2020")
        row["effect"]["seed_output_omnibus_p"] = 0.01
        with self.assertRaisesRegex(ValueError, "Source overall seed test changed"):
            synthesize(tmp)


if __name__ == "__main__":
    unittest.main()
