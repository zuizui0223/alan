#!/usr/bin/env python3
"""Dryad moth source grain maths and frozen corrected source file identity."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_moth23_Dryad_corrected_source_grain_20261114 import (
 FILES,median,BASE
)
class Moth23Original(unittest.TestCase):
 def test_pinned_corrected_data(self):
  self.assertEqual(FILES["Activity_MainData.csv"][0],"4507228")
  self.assertEqual(FILES["Activity_MainData.csv"][1],77107372)
  self.assertEqual(FILES["Activity_MainData.csv"][2],
   "e08798c04d5708e5ea51e0ceda8f4ea64f192b38656726c1f7b9a94305f18c04")
 def test_no_frame_independence_presumed(self):
  self.assertEqual(median([4,2,3]),3)
  self.assertEqual(median([4,2,3,5]),3.5)
  self.assertTrue(BASE.startswith("https://datadryad.org/api/v2/files/"))
if __name__=="__main__":unittest.main()
