"""Recompute the post hoc n=3 Sullivan contrast from its original-source frozen blocks."""
import copy
import csv
import importlib.util
import json
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "database/alan_experiments/arthropod_expansion_20261122/sullivan2019_source_reanalysis"
SCRIPT = STUDY / "colour_balanced_block_test_20261124.py"
SOURCE = STUDY / "dryad_provenance_and_frozen_block_results.json"
PREVIOUS = STUDY / "within_experiment_cross_endpoint_ratios.csv"

spec = importlib.util.spec_from_file_location("sullivan_colour_balanced", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestColourBalancedSullivan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(SOURCE.read_text(encoding="utf-8"))
        cls.out = module.analyze(cls.data)

    def test_three_physical_blocks_not_four_colours_as_studies(self):
        for key in ("TerrestrialArthropod", "WebCount"):
            result = self.out["comparisons"][key]
            self.assertEqual(result["n_independent_habitat_blocks"], 3)
            self.assertEqual(len(result["block_ratios"]), 3)
            self.assertTrue(all(x > 1 for x in result["block_ratios"]))
            self.assertAlmostEqual(result["block_signflip_two_sided_p"], .25)
            self.assertEqual(len(result["leave_one_colour_out"]), 4)
            self.assertTrue(all(v["colour_balanced_ratio_of_RRs"] > 1
                                for v in result["leave_one_colour_out"].values()))

    def test_frozen_omnibus_numbers_and_multiplicity(self):
        a = self.out["comparisons"]["TerrestrialArthropod"]
        w = self.out["comparisons"]["WebCount"]
        self.assertAlmostEqual(a["colour_balanced_ratio_of_RRs"], 2.4739657387985243, places=10)
        self.assertAlmostEqual(a["approx_t95_ratio_interval"][0], 0.3817194652222584, places=10)
        self.assertAlmostEqual(a["approx_t95_ratio_interval"][1], 16.034043412444863, places=10)
        self.assertAlmostEqual(w["colour_balanced_ratio_of_RRs"], 7.322253300306327, places=10)
        self.assertAlmostEqual(w["approx_t95_ratio_interval"][0], 1.4654010422979726, places=10)
        self.assertAlmostEqual(w["approx_t95_ratio_interval"][1], 36.58752235481543, places=10)
        self.assertAlmostEqual(w["approx_t_two_sided_p"], 0.03350834516975072, places=10)
        self.assertAlmostEqual(w["exploratory_two_contrast_Holm_p"], 0.0670166903395014, places=10)

    def test_each_colour_matches_previously_published_derived_csv(self):
        ep = {e: module.endpoint_blocks(self.data, e) for e in module.ENDPOINTS}
        with PREVIOUS.open(newline="", encoding="utf-8") as h:
            rows = list(csv.DictReader(h))
        self.assertEqual(len(rows), 8)
        for row in rows:
            result = module.block_comparison(ep["AquaticEmergent"],
                                             ep[row["comparison_response"]],
                                             (row["lamp_treatment"],))
            self.assertTrue(math.isclose(result["average_ln_ratio_of_RRs"],
                                         float(row["difference_lnRR"]), abs_tol=1e-12))
            self.assertTrue(math.isclose(result["se_across_blocks"],
                                         float(row["SE_block_paired"]), abs_tol=1e-12))

    def test_rejects_corrupted_source_summary(self):
        corrupted = copy.deepcopy(self.data)
        corrupted["provisional_block_contrasts"]["results"]["AquaticEmergent"]["estimands"]["LEDI"]["lnRR"] += .1
        with self.assertRaises(ValueError):
            module.analyze(corrupted)

    def test_deterministic(self):
        self.assertEqual(self.out, module.analyze(self.data))


if __name__ == "__main__":
    unittest.main()
