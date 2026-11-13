#!/usr/bin/env python3
"""Only genuine independent positive GloBI flower interactions count as evidence.
No source-missing field study plant can be classified as biologically day only.
"""
import unittest
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_genuine_diel_GloBI_v011_flower_guild_overlap_20261113 import (
    FLAGS_SHA,SOURCE_SHA,TIERS,compare_positive_guild_documentation,
    MIN_GROUP_STUDIES,MIN_GROUP_SPECIES,LIGHT
)


class GloBITrueFieldDielSourceTests(unittest.TestCase):
    def test_original_scanner_source_pins_and_methods(self):
        self.assertEqual(len(SOURCE_SHA),64)
        self.assertEqual(len(FLAGS_SHA),64)
        self.assertEqual(MIN_GROUP_STUDIES,25)
        self.assertEqual(MIN_GROUP_SPECIES,20)
        self.assertIn("non_inat_visits_source",TIERS)

    def test_missing_butterfly_record_means_unknown_not_biological_absence(self):
        f=pd.DataFrame({
           "study_key":["a","b","c"],"species":["A b","C d","E f"],
           "night_share":[.85,.35,.6],"plant_family":["f1","f2","f3"],
           "bloom_period_simple":["night","day","unknown"],
           LIGHT:[0,.4,1.3],
           "moth_record":[1,0,0],"butterfly_record":[0,1,0]
        })
        r=compare_positive_guild_documentation(f,"synthetic_positive_only",
                       "moth_record","butterfly_record")
        self.assertEqual(
           r["observed_positive_guild_documentation_counts"]
            ["moth_positive_butterfly_NOT_documented"]["distinct_original_plant_species"],1)
        self.assertEqual(
           r["can_infer_real_moth_vs_butterfly_original_reproductive_contrast"]["status"],
           "NOT_ESTIMABLE_INDEPENDENT_SOURCE_PLANT_AND_PAPER_GUILD_SUPPORT")
        self.assertNotIn("p",r["can_infer_real_moth_vs_butterfly_original_reproductive_contrast"])
        self.assertTrue(r["biological_limits"]["GloBI_butterfly_NOT_documented_does_NOT_mean_no_daytime_pollinators"])

if __name__=="__main__":
    unittest.main()
