#!/usr/bin/env python3
"""Synthetic no-p-hacking unit controls for true original study reporting bias."""
import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original301_observer_OA_sampling_bias_20261114 import (
    difference,require,SOURCE_SHA,AUTHOR_GITBLOB
)


class Original301ObserverSelection(unittest.TestCase):
    def test_study_level_descriptive_no_individual_pseudo_replication(self):
        a=np.array([.8]*30,dtype=float)
        b=np.array([.2]*25,dtype=float)
        x=difference(a,b)
        self.assertEqual(x["original_independent_studies_groupA"],30)
        self.assertEqual(x["original_independent_studies_groupB"],25)
        self.assertAlmostEqual(x["groupA_minus_B_mean_study_median_difference"],.6)
        self.assertAlmostEqual(x["descriptive_study_bootstrap_95CI"][0],.6)
        self.assertNotIn("p_value",x)
        self.assertEqual(x["study_block_bootstrap_draws"],3000)

    def test_not_enough_independent_papers_fails_closed(self):
        with self.assertRaises(ValueError):
            difference([.2]*21,[.9]*8)

    def test_original_source_gates_are_pinned(self):
        self.assertEqual(len(SOURCE_SHA),64)
        self.assertEqual(len(AUTHOR_GITBLOB),40)
        with self.assertRaises(ValueError):
            require(False,"OA metadata cannot stand for effective moth pollen")


if __name__=="__main__":
    unittest.main()
