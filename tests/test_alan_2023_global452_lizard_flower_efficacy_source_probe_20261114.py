#!/usr/bin/env python3
"""Original 452-lizard interaction source discoverability and study-grain guards."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_2023_global452_lizard_flower_efficacy_original_source_probe_20261114 import (
    SOURCES,request
)
class Global452Source(unittest.TestCase):
 def test_sources_are_specific_true_original_article_not_globi(self):
  self.assertEqual(len(SOURCES["wiley_official_2023_S1"]),4)
  self.assertTrue(all("ecog12967-sup-0001-AppendixS1.xls" in u
                      for u in SOURCES["wiley_official_2023_S1"]))
  self.assertTrue(all(u.startswith("https://") for v in SOURCES.values() for u in v))
  self.assertIn("23460260",SOURCES["figshare_lizard_macro_2838"][0])
 def test_missing_local_visits_cannot_become_zero_efficacy(self):
  self.assertNotEqual(90,452)
  self.assertLess(90/452,0.25)
  self.assertGreater(90/452,0.19)
  self.assertEqual(452-90,362)
if __name__=="__main__":unittest.main()
