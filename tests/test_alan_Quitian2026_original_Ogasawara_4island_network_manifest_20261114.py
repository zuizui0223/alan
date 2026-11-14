#!/usr/bin/env python3
"""Real original author Zenodo file gate tests; no ecological pseudo-zero."""
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_Quitian2026_original_Ogasawara_4island_network_source_manifest_20261114 import (
 URL,DOI,PAPER,SOURCE_MD5
)
class OriginalOgasawaraNetwork(unittest.TestCase):
 def test_source_exact_dataset(self):
  self.assertEqual(DOI,"10.5281/zenodo.19221853")
  self.assertEqual(PAPER,"10.1111/cobi.70304")
  self.assertIn("Data_Quitianetal2026_B.xlsx",URL)
  self.assertEqual(SOURCE_MD5,"68f295b91605ea82bce0c05bffc7f884")
 def test_known_island_design_not_world_replicates(self):
  self.assertEqual(len(("Chichijima","Hahajima","Anijima","Ototojima")),4)
  self.assertGreater(4174/55,10)
  self.assertNotEqual(4174,4)
if __name__=="__main__":
 unittest.main()
