import json
from pathlib import Path
import unittest

REG = json.loads((Path(__file__).resolve().parents[1] / "config/alan_A_reproductive_screen_20261115.json").read_text(encoding="utf-8"))
KNOP = next(x for x in REG["entries"] if x["id"]=="knop2017")["source_audit"]

class TestKnopRealSourceDesign(unittest.TestCase):
    def test_light_assigned_fruit_vs_visits_not_same(self):
        self.assertEqual((KNOP["night_network_lit_sites"],KNOP["night_network_dark_sites"]),(7,7))
        self.assertEqual((KNOP["fruit_lit_sites"],KNOP["fruit_dark_sites"]),(5,5))
        self.assertEqual(KNOP["fruit_field_sites_total"],10)
        self.assertEqual(KNOP["matched_site_pairs"],5)
    def test_true_nested_plants_and_heads(self):
        self.assertEqual(KNOP["total_fruit_plants"],100)
        self.assertEqual(KNOP["plants_per_lit_x_bagging_cell"],25)
        self.assertEqual(KNOP["capitula_per_plant"],2)
        self.assertEqual(KNOP["total_measured_capitula"],200)
        self.assertEqual(KNOP["capitula_per_lit_x_bagging_cell"],50)
        self.assertEqual(4*KNOP["capitula_per_lit_x_bagging_cell"],KNOP["total_measured_capitula"])
        self.assertEqual(KNOP["total_fruit_plants"],10*KNOP["plants_per_site"])
    def test_not_diel_dependency_and_no_se(self):
        self.assertFalse(KNOP["full_day_night_exclusion"])
        self.assertTrue(KNOP["bagging_does_not_separate_diel_guilds"])
        self.assertFalse(KNOP["original_model_variance_available_to_meta"])
        self.assertFalse(KNOP["full_interaction_numeric_beta_se_recovered"])
    def test_author_p_values_are_not_interaction(self):
        self.assertAlmostEqual(KNOP["light_negative_open_fruits_p"],.014)
        self.assertAlmostEqual(KNOP["light_effect_on_bagged_fruits_p"],.140)
        self.assertFalse(next(x for x in REG["entries"] if x["id"]=="knop2017")["verified_joint_effect_moderator_se"])

if __name__ == "__main__":
    unittest.main()
