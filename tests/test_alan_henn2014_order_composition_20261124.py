"""Source-derived six-order compositional check for Henn 2014; n = five streams."""
import copy
import csv
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124"
spec=importlib.util.spec_from_file_location("henn2014_composition", SOURCE/"henn2014_order_composition.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestHennOrderComposition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (SOURCE/"henn2014_derived_order_by_stream.csv").open(encoding="utf-8",newline="") as fh:
            cls.rows=list(csv.DictReader(fh))
        cls.result=module.analyze(cls.rows)

    def test_five_genuine_streams_and_source_mismatch_are_kept(self):
        a=self.result
        self.assertEqual(a["n_independent_stream_pairs"],5)
        self.assertEqual(len(a["stream_results"]),5)
        self.assertEqual(a["sum_source_family_group_counts"],{"lit":1231,"dark":1961,"total":3192})
        self.assertEqual(a["sum_original_table_total_row"],{"lit":1229,"dark":1961,"total":3190})
        self.assertEqual(a["family_minus_table_total"],{"lit":2,"dark":0})
        self.assertEqual(a["meta_programme_increment"],0)

    def test_compositional_ratio_exact(self):
        a=self.result
        self.assertAlmostEqual(a["equal_stream_compositional_ratio_of_RRs"],1.2846373445392218,places=12)
        self.assertAlmostEqual(a["descriptive_t4_95_interval"][0],.8940489154931339,places=8)
        self.assertAlmostEqual(a["descriptive_t4_95_interval"][1],1.8458644469967567,places=8)
        self.assertEqual(a["posthoc_two_sided_signflip_p"],.0625)
        self.assertTrue(a["all_five_ratios_above_one"])
        self.assertAlmostEqual(a["stream_results"][0]["composition_ratio_of_RRs"],1.58974358974359,places=12)
        self.assertAlmostEqual(a["stream_results"][3]["composition_ratio_of_RRs"],1.001082251082251,places=12)

    def test_source_mismatch_bound_and_leave_one_out(self):
        a=self.result
        lo,hi=a["hypothetical_san_marcos_two_count_allocation_bound"]
        self.assertAlmostEqual(lo,1.283718432765531,places=12)
        self.assertAlmostEqual(hi,1.286854160740724,places=12)
        vals=list(a["leave_one_stream_out"].values())
        self.assertTrue(all(1.16<v<1.38 for v in vals))
        q=a["four_abundant_families"]
        self.assertEqual((q["lit"],q["dark"]),(909,1459))
        self.assertAlmostEqual(q["share_of_source_aggregate_count_deficit"],550/730,places=12)
        self.assertAlmostEqual(q["share_of_all_captured_individuals"],2368/3192,places=12)

    def test_invalid_family_aggregation_stops_instead_of_adjusting(self):
        rows=copy.deepcopy(self.rows)
        rows[0]["SanMarcos_lit"]="90"
        with self.assertRaises(ValueError):
            module.analyze(rows)
        rows=copy.deepcopy(self.rows)
        rows[1]["order"]=rows[0]["order"]
        with self.assertRaises(ValueError):
            module.analyze(rows)

    def test_no_pseudocounts_and_same_published_totals(self):
        self.assertEqual(self.result,module.analyze(self.rows))
        with self.assertRaises(ValueError):
            module.one_river(0,1,2,3)


if __name__=="__main__":
    unittest.main()
