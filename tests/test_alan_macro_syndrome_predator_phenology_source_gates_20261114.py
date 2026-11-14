#!/usr/bin/env python3
"""Tests on actual versioned public source results: no fake meta-effect."""
import tempfile
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_macro_syndrome_predator_phenology_source_gates_20261114 import (
    finite_ci, require, main
)

ROOT=Path(__file__).resolve().parents[1]


class MacroNightPredatorPhenologyGates(unittest.TestCase):
    def test_real_original_signed_aggregate_source_gate(self):
        with tempfile.TemporaryDirectory() as d:
            result=main(ROOT,Path(d)/"source_audit.json")
        a=result["original_source_aggregate_replay"]
        self.assertEqual(a["island_macro_observed_photo_opportunity"]["independent_islands"],562)
        self.assertEqual(a["island_mainland_photo_alanspatial"]["island_siteyears"],933)
        self.assertEqual(a["island_mainland_photo_alanspatial"]["mainland_siteyears"],9765)
        self.assertEqual(a["animal_photo_predator_and_moth"]["matched_same_island_cell_years"],15255)
        self.assertEqual(a["reproductive_field_night_only_share_vs_historical_light"]["pairs"],301)
        self.assertEqual(len(result["conditional_hypotheses_identification"]),4)
        self.assertEqual([x["id"] for x in result["conditional_hypotheses_identification"]],
                         ["H1","H2","H3","H4"])
        self.assertEqual(result["conditional_hypotheses_identification"][2]["direct_complete_local_flower_pollinator_phenophase_pairs_under_light"],0)
        self.assertTrue(result["safe_inference"]["local_phenology_and_syndrome_shift_still_missing_for_macro_4way"])
        self.assertNotIn("global_joint_interaction_p",result)

    def test_cross_source_false_PL64_p_transcription_is_explicit_not_silently_copied(self):
        with tempfile.TemporaryDirectory() as d:
            result=main(ROOT,Path(d)/"source.json")
        x=result["old_source_record_provenance_bug"]
        self.assertAlmostEqual(x["actually_readme_PL64_p_from_direct_original_replay"],.043993409104236464,places=10)
        self.assertAlmostEqual(x["incorrectly_copied_v011_PL30_p"],.4799746115289204,places=10)
        self.assertNotAlmostEqual(x["actually_readme_PL64_p_from_direct_original_replay"],
                                 x["value_currently_written"],places=4)
        self.assertTrue(x["prior_original_file_not_overwritten"])

    def test_invalid_confidence_bounds_not_citable(self):
        with self.assertRaises(ValueError):
            finite_ci({"ci95":[1,-1]},"ci95")
        with self.assertRaises(ValueError):
            finite_ci({"ci95":[float("nan"),1]},"ci95")
        self.assertEqual(finite_ci({"ci95":[-1,1]},"ci95"),[-1,1])


if __name__=="__main__":
    unittest.main()
