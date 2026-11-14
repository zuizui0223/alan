#!/usr/bin/env python3
"""Real source network join, effective number of source visits and no false biology."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_Quitian2026_Ogasawara4island_visitation_functional_portfolio_20261114 import (
 entropy_effective,api,COLUMNS,ISLANDS,SEASONS,model
)
class OgasawaraSourceMath(unittest.TestCase):
 def test_partner_number_is_visit_entropy_not_effective_pollen(self):
  self.assertAlmostEqual(entropy_effective({"A":5,"B":5}),2,places=10)
  self.assertAlmostEqual(entropy_effective({"A":5,"B":0}),1,places=10)
  self.assertIsNone(entropy_effective({}))
 def test_honeybee_string_and_other_insect_not_native(self):
  self.assertTrue(api("Apis mellifera"))
  self.assertTrue(api(" A. mellifera "))
  self.assertFalse(api("Ceratina boninensis"))
  self.assertEqual(len(COLUMNS),11)
  self.assertEqual(len(ISLANDS),4)
  self.assertEqual(len(SEASONS),3)
if __name__=="__main__":unittest.main()
