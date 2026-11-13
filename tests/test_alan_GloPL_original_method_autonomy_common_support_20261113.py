#!/usr/bin/env python3
"""GloPL method selfing x historical ALAN must never be called same-paper effect
without genuine matched original treatment support."""
import sys
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_GloPL_original_method_autonomy_common_support_20261113 import (
    within_research_light_support,sup_only_selfing_direct,
)

class RealSourceStudyMethodAuditTests(unittest.TestCase):
    def test_two_methods_in_different_papers_not_matched_design(self):
        q=pd.DataFrame({
          "study_key":["p1","p1","p2","p2"],
          "autonomy":[0,1,0,1],
          "PL_Effect_Size_Type2":["Sup","Sup","Bagout","Bagout"],
          "pcnl_prior_mean_log1p":[0.1,.15,1.1,1.3],
          "original_row_id":[1,2,3,4]})
        result=within_research_light_support(q)
        self.assertEqual(result["source_publications_with_both_selfing_types"],2)
        self.assertEqual(result["publications_with_both_methods_both_selfing"],0)
        self.assertEqual(result["source_publications_both_selfing_types_AND_light_range_at_least_003"],2)

    def test_no_unsupported_pseudo_effect_for_rare_sup_source(self):
        q=pd.DataFrame({
          "PL_Effect_Size_Type2":["Sup"]*4,
          "study_key":["a","b","c","d"],
          "autonomy":[0,0,1,1],
          "PL_log_ratio":[.1,.2,.3,.4],
          "pcnl_prior_mean_log1p":[0.,.1,.2,.3]
        })
        m=sup_only_selfing_direct(q)
        self.assertEqual(m["status"],"NOT_ESTIMABLE_SUP_STUDY_SOURCE_SUPPORT")
        self.assertNotIn("beta_PL_log_reproductive_supplementation_ratio_for_autonomous_vs_absent",m)

if __name__=="__main__":
    unittest.main()
