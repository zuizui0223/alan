import importlib.util
import json
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_B_original_DvN_stage_contrasts_20261115.py"
spec=importlib.util.spec_from_file_location("stages",P)
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)

class TestStageContrast(unittest.TestCase):
    def setUp(self):
        self.rows=s.parse_rows()
    def test_42_original_rows(self):
        self.assertEqual(len(self.rows),42)
    def test_14_stage_pairs_and_strict_10(self):
        x=s.match_and_gate(self.rows)
        self.assertEqual(x["all_stage_pairs"],14)
        self.assertEqual(x["strict_stage_pairs"],10)
        self.assertEqual(x["strict_original_studies"],6)
        self.assertEqual(x["rejected_pair_count"],4)
    def test_stage_differences_descriptive_not_meta(self):
        x=s.match_and_gate(self.rows)
        self.assertEqual(x["strict_seed_marginal_larger_count"],10)
        self.assertEqual(x["all_pair_seed_marginal_larger_count"],14)
        self.assertIsNone(x["global_meta_beta_V"])
        self.assertFalse(x["same_initial_flower_total_seed_fitness_identified"])
    def test_mismatched_french_original_quality_excluded(self):
        x=s.match_and_gate(self.rows)
        reason=next(r for r in x["rejected_reasons"] if r["study_ID"]=="SmallFruit_2025")
        self.assertIn("MISMATCH",reason["reason"])
        self.assertFalse(any(p["study_ID"]=="SmallFruit_2025" for p in x["strict"]))
    def test_morse_4_strata_one_publication(self):
        x=s.match_and_gate(self.rows)
        m=[a for a in x["strict"] if a["study_ID"]=="Morse_1983"]
        self.assertEqual(len(m),4)
        self.assertTrue(all(a["V_seed_minus_V_fruit"]>0 for a in m))
    def test_reject_changed_reference_day_open(self):
        rows=[dict(r) for r in self.rows]
        seed=next(r for r in rows if r["study_ID"]=="Morse_1983" and r["treatment_effectiveness_metric"]=="seed set")
        seed["day_value"]="1234"
        with self.assertRaises(AssertionError):s.match_and_gate(rows)

if __name__=="__main__":
    unittest.main()
