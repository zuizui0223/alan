import csv
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124"
spec=importlib.util.spec_from_file_location("henn2014_stream_counts", SOURCE/"henn2014_table1_analysis.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestHenn2014(unittest.TestCase):
    def test_original_table_recovery(self):
        with (SOURCE/"henn2014_table1_five_stream_counts.csv").open(newline="",encoding="utf-8") as fh:
            result=module.aggregate(list(csv.DictReader(fh)))
        self.assertEqual(result["table1_count_sums"], {"lit":1229,"dark":1961,"combined":3190})
        self.assertAlmostEqual(result["grand_total_lit_over_dark"], .6267210606833248,places=12)
        self.assertAlmostEqual(result["stream_equal_geometric_mean_RR"],.9998109779357996,places=12)
        self.assertAlmostEqual(result["descriptive_t4_interval_stream_RR"][0],.22629854340747577,places=12)
        self.assertAlmostEqual(result["descriptive_t4_interval_stream_RR"][1],4.41727099321629,places=12)
        self.assertEqual(result["unadjusted_32_signflip_two_sided_p_sensitivity"],1.0)
        self.assertAlmostEqual(result["leave_one_stream_out_equal_weight_ratio"]["Barton Creek"],.591838444225893,places=12)
        self.assertAlmostEqual(result["size_strata_exploratory"]["large"]["geometric_mean_stream_RR"], .5610071966833595,places=12)
        self.assertEqual(result["meta_adopted"],False)

    def test_rejects_corruption(self):
        with (SOURCE/"henn2014_table1_five_stream_counts.csv").open(newline="",encoding="utf-8") as fh:
            rows=list(csv.DictReader(fh))
        rows[0]["lit_count"]="58"
        with self.assertRaises(ValueError):
            module.aggregate(rows)
        rows[0]["lit_count"]="57"
        rows[1]["stream"]=rows[0]["stream"]
        with self.assertRaises(ValueError):
            module.aggregate(rows)


if __name__=="__main__":
    unittest.main()
