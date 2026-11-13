#!/usr/bin/env python3
"""Source-safe numerical guardrails for genuine day/night study vs island PCNL."""
import sys
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import (
    source_regression, signed_contrast, blob_hash, canon, KEY,
)

class SourceControls(unittest.TestCase):
    def test_missing_names_are_not_positive(self):
        x=canon(pd.Series([" Salix   caprea ","","Some sp"]))
        self.assertEqual(x.tolist(),["Salix caprea","","Some sp"])

    def test_rare_night_response_returns_unestimated_not_pseudo_p(self):
        sites=pd.DataFrame({
          "island":["A","A","B"],"latbin":[1,2,1],"lonbin":[1,1,3],
          "year":[2023,2023,2023],"all_accepted_recorded_species":[3,3,4]})
        plants=pd.DataFrame({
          "island":["A","A","B"],"latbin":[1,2,1],"lonbin":[1,1,3],
          "year":[2023,2023,2023],"night_stronger":[True,False,False],
          "day_stronger":[False,True,False]
        })
        m=source_regression(plants,sites,"night_stronger")
        self.assertEqual(m["status"],"NOT_ESTIMABLE_INDEPENDENT_ISLAND_SUPPORT")
        self.assertNotIn("fit",m)
        v=signed_contrast(plants,sites,"night_stronger","day_stronger")
        self.assertEqual(v["status"],"NOT_ESTIMABLE_SIGNED_CONTRAST_INDEPENDENT_ISLAND_SUPPORT")
        self.assertNotIn("model",v)

    def test_nonmatching_original_blob_refuses(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/"x.csv"
            path.write_text("night,value\nmock,1\n")
            self.assertNotEqual(blob_hash(path),"39a82d966afcdde670603add1c881ace8791ffe9")

if __name__=="__main__":
    unittest.main()
