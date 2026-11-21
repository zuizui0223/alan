"""Fail-closed checks for source-fixed ALAN wet-habitat candidate-mechanism triage."""
import math
import unittest
from scripts.analysis import alan_wet_candidate_triage_20261117 as triage


class TestWetCandidateTriage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = triage.inspect()

    def test_original_primary_effect_not_reselected(self):
        r = self.r["primary_existing_moisture_difference"]
        self.assertEqual((r["n_candidate"], r["n_reference"], r["regions"]),
                         (158, 1419, 12))
        self.assertTrue(math.isclose(r["unadjusted"], 0.5922104455502863))
        self.assertTrue(math.isclose(r["adjusted"], 0.5338399813053877))
        self.assertGreater(r["adjusted_interval"][0], 0)

    def test_heterogeneous_regions_and_broad_modifiers_stay_negative(self):
        region = self.r["region_transport_diagnostic"]
        self.assertEqual((region["positive"], region["negative"]), (7, 5))
        self.assertTrue(math.isclose(region["top_two_regions_candidate_share"], 94/158))
        m = self.r["previous_macro_falsifiers"]
        self.assertEqual(m["extra_trait_count"], 34)
        self.assertGreater(m["best_adjusted_q_extra_traits"], .05)
        self.assertLess(m["regional_AUC_plus_34_traits"], m["regional_AUC_baseline"])
        self.assertEqual(m["site_modifiers_count"], 12)
        self.assertGreater(m["best_site_modifier_q"], .05)

    def test_favorable_four_taxon_sign_does_not_promote_a_meta(self):
        s = self.r["swiss_same_experiment_moisture_pilot"]
        self.assertEqual((s["experiment_count"], s["all_plant_taxa"],
                          s["wet_indicator_known_species"]), (1, 21, 4))
        self.assertTrue(math.isclose(s["spearman_r_descriptive"], -1.0))
        self.assertTrue(math.isclose(s["spearman_exact_permutation_two_sided"], 2/24))
        self.assertEqual(s["drier_species_with_at_least_two_pairs"], 0)
        self.assertFalse(s["eligible_as_cross_species_effect_meta"])
        self.assertFalse(self.r["meta_readiness"]["homogeneous_multi_study_lit_dark_plant_response_with_wet_dry_moderator"])


if __name__ == "__main__":
    unittest.main()
