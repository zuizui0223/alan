import importlib.util
import math
from pathlib import Path
import unittest

SRC = Path(__file__).resolve().parents[1] / "scripts/analysis/alan_A_effect_type_gate_20261115.py"
SPEC = importlib.util.spec_from_file_location("alan_a_effect", SRC)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

class TestALANEffectTypeGate(unittest.TestCase):
    def test_li_plant_binary_not_flower(self):
        e = MODULE.publication_compatible_examples()["studies"][0]
        self.assertAlmostEqual(e["light_dark_ratio"],(20/27)/(8/26))
        self.assertAlmostEqual(e["ln_RR"],math.log((20/27)/(8/26)))
        self.assertFalse(e["primary_initial_flower_endpoint"])
        self.assertIsNone(e["ln_RR_SE_independent_clustering"])
    def test_lu_no_guessed_se(self):
        e=MODULE.publication_compatible_examples()["studies"][1]
        self.assertAlmostEqual(e["ln_ratio_of_author_percent"],math.log(.44))
        self.assertIsNone(e["ln_RR_SE"])
    def test_giavi_or_not_rr(self):
        e=MODULE.publication_compatible_examples()["studies"][2]
        self.assertAlmostEqual(e["derived_odds_ratio"],math.exp(-1.37))
        self.assertIsNone(e["ln_RR"])
    def test_boom_no_false_significance(self):
        e=MODULE.publication_compatible_examples()["studies"][3]
        self.assertAlmostEqual(e["ln_RR_approx_product"],math.log(.5150))
        self.assertIsNone(e["ln_RR_SE"])
    def test_reject_zero_controls(self):
        with self.assertRaises(ValueError):
            MODULE.log_response_ratio(1.0,0.0)
    def test_no_beta(self):
        self.assertIsNone(MODULE.publication_compatible_examples()["confirmed_A_beta_V"])

if __name__=="__main__":
    unittest.main()
