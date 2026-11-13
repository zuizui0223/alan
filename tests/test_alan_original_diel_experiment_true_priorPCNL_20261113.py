#!/usr/bin/env python3
"""Fail-closed study-level pairing and original prior light chronology."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_diel_experiment_true_priorPCNL_20261113 import (
    YEARS,SHARDS,PAIR_N,STUDY_N,SPECIES_N,AUTHOR_BLOB,
    shard_years,wls_study
)


class IndependentDayNightExperimentSourceTests(unittest.TestCase):
    def test_satellite_year_partition_exactly_once(self):
        all_years=[x for i in range(SHARDS) for x in shard_years(i)]
        self.assertEqual(sorted(all_years),list(range(1992,2022)))
        self.assertEqual(len(all_years),len(set(all_years)))
        self.assertEqual(len(all_years),30)

    def test_old_outcome_source_lock(self):
        self.assertEqual((PAIR_N,STUDY_N,SPECIES_N),(146,87,94))
        self.assertEqual(len(AUTHOR_BLOB),40)

    def test_no_bogus_p_from_few_independent_studies(self):
        z=pd.DataFrame({
           "study":["one"]*10,
           "species":["Plant sp"]*10,
           "geo10":["0_0"]*10,
           "pcnl_prior_mean_log1p":np.linspace(.1,1.,10),
           "nightshare_pp":np.linspace(10,80,10)})
        r=wls_study(z,"nightshare_pp ~ pcnl_prior_mean_log1p")
        self.assertEqual(r["status"],"NOT_ESTIMABLE_STUDY_LIGHT_SUPPORT")
        self.assertNotIn("beta_nightshare_pp_per_mean_logPCNL",r)

    def test_study_sample_zero_no_p_inferred(self):
        z=pd.DataFrame({
          "study":["a"]*3+["b"]*3,"species":["A"]*6,
          "geo10":["10_10"]*6,
          "pcnl_prior_mean_log1p":[.1,.2,.3,.1,.2,.3],
          "nightshare_pp":[45,55,60,44,50,60]})
        r=wls_study(z,"nightshare_pp ~ pcnl_prior_mean_log1p")
        self.assertEqual(r["status"],"NOT_ESTIMABLE_STUDY_LIGHT_SUPPORT")


if __name__=="__main__":
    unittest.main()
