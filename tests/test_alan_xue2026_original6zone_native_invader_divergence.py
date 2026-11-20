"""Biological source/physical-zone math, no 77-cylinder pseudo-replication."""
import math, unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_xue2026_original6zone_native_invader_divergence import (
    read_source, validate_zones, evaluate_labels, species_effects,
    community_description, git_blob, SRC
)

class TestOriginalXueZones(unittest.TestCase):
    def test_pinned_original_source_grid_exactly_six_zones(self):
        x=read_source()
        p,l=validate_zones(x["zone"])
        self.assertEqual(p.shape,(6,4))
        self.assertEqual(int(x["zone"].n_surviving_cylinders.sum()),77)
        self.assertEqual(int((l=="ALAN").sum()),3)
        self.assertEqual(int((l=="No ALAN").sum()),3)

    def test_actual_source_index_cannot_claim_77_independent_lamps(self):
        x=read_source()
        p,l=validate_zones(x["zone"])
        y=evaluate_labels(p,l,reps=149)
        self.assertEqual(y["independent_lit_light_zones"],3)
        self.assertEqual(y["independent_dark_light_zones"],3)
        self.assertEqual(y["exact_exchangeable_three_of_six_zone_assignments"],20)
        self.assertAlmostEqual(y["exact_zone_label_two_sided_p_if_exchangeable"],.1)
        self.assertGreater(y["invasive_to_native_relative_response_ratio"],1.6)
        self.assertTrue(y["all_three_lit_zone_indices_exceed_every_dark_zone_index"])
        self.assertTrue(y["not_proven_original_randomization_p"])

    def test_shared_light_zone_species_covariance_is_preserved(self):
        x=read_source()
        p,l=validate_zones(x["zone"])
        z=species_effects(p,l,x["effects"],x["covariance"])
        self.assertAlmostEqual(
           z["source_invasive_lnRR_minus_native_average_lnRR"],.4652787951402042,places=8)
        self.assertGreater(
           z["within_source_lighting_zone_covariance_delta_approx_variance"],0.)

    def test_zone_duplicate_cannot_be_added_as_replication(self):
        x=read_source()["zone"]
        fake=pd.concat([x,x.iloc[[0]]],ignore_index=True)
        with self.assertRaisesRegex(RuntimeError,"24 unique"):
            validate_zones(fake)

    def test_original_community_has_missing_fixture_assignment_stop(self):
        z=community_description(read_source()["community"])
        self.assertTrue(z["underlying_real_independent_lamp_count_UNRESOLVED"])
        self.assertAlmostEqual(
          z["summary_by_original_light_arm"]["ALAN"]["native_suppression_relative_to_same_light_no_snail"],
          .87225,places=3)
        self.assertFalse(z.get("causal_grazer_mediation_confirmed",False))

if __name__=="__main__":
    unittest.main()
