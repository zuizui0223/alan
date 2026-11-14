#!/usr/bin/env python3
"""Two true original field regions: mathematics and observational-inference guards."""
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_two_region_original_guild_flowercover_20261114 import (
    jsbits,lograte,summarize,ORDERS
)

class FunctionalDayGuildSource(unittest.TestCase):
    def test_Jensen_Shannon_original_three_order_axis(self):
        self.assertAlmostEqual(jsbits([10,20,30],[10,20,30]),0,places=12)
        self.assertAlmostEqual(jsbits([1,0,0],[0,1,0]),1,places=12)
        self.assertAlmostEqual(jsbits([0,2,1],[0,2,1]),0,places=12)
        with self.assertRaisesRegex(ValueError,"No biological"):
            jsbits([0,0,0],[1,0,0])
    def test_exposure_adjustment_not_faked_as_pollen(self):
        self.assertAlmostEqual(lograte(20,10,2,1),0,places=8)
        self.assertAlmostEqual(lograte(20,10,1,1),__import__("math").log(20.5/10.5),places=12)
    def test_true_field_block_aggregate_structure(self):
        records=[]
        for i in range(1,6):
            counts={}
            for k in ORDERS:
                counts[k]={
                   "log_lit_dark_per_survey":lograte(10+i,10,1,1),
                   "log_lit_dark_per_recorded_flower_cover":lograte(10+i,10,2,1),
                   "lit_original_visits":10+i,
                   "dark_original_visits":10}
            records.append({"guilds":counts,"three_order_JS_divergence_bits":0.15})
        a=summarize(records)
        self.assertEqual(a["Diptera"]["independent_sites_or_blocks"],5)
        self.assertEqual(a["Diptera"]["original_sites_where_lit_visit_rate_higher"],5)
        self.assertTrue(a["Diptera"]["no_claimed_causal_pollen_service_or_global_meta_effect"])
        self.assertAlmostEqual(a["three_order_JS_composition"]["median_JS_divergence_bits"],.15)
        b=summarize(records,cover=True)
        self.assertLess(b["Diptera"]["equal_site_geometric_mean_lit_dark_VISIT_RATE_ratio"],
                        a["Diptera"]["equal_site_geometric_mean_lit_dark_VISIT_RATE_ratio"])

if __name__=="__main__":
    unittest.main()
