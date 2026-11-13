#!/usr/bin/env python3
"""Test study-level light common-support and refusal of fake source traits."""
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_DvN_real_field_light_study_influence_bloom_20261113 import (
   actual_comparable_light, original_trait_bloom, BASELINE,
)


class RealPairedFieldSupportTests(unittest.TestCase):
    def test_source_geography_study_does_not_invent_independent_experiment(self):
        d=pd.DataFrame({
            "study_key":["p1","p1","p2","p2"],
            "source_pair_id":[1,2,3,4],
            "species":["spA","spA","spB","spB"],
            "lat":[10,10,10,10],"lon":[70,70,70,70],
            "pcnl_prior_mean_log1p":[0.,.09,.0,.3]
        })
        z=actual_comparable_light(d)
        self.assertEqual(z["source_light_zero_pairs"],2)
        self.assertEqual(z["source_light_positive_pairs"],2)
        self.assertEqual(z["within_real_published_study"]
                         ["studies_with_both_zero_and_positive_prior_satellite_PCNL"],2)
        self.assertEqual(z["10degree"]["regions_containing_both_zero_and_positive_prior_PCNL"],1)
        self.assertEqual(z["10degree"]["regions_both_arms_at_least_2_independent_publications"],1)

    def test_reject_unverified_original_botanical_source(self):
        with TemporaryDirectory() as d:
            p=Path(d)/"wrong.csv"
            p.write_text("species,bloom_period_simple\nfake,night\n")
            with self.assertRaises(ValueError):original_trait_bloom(p)

    def test_real_paired_field_estimator_benchmark_preserved(self):
        self.assertAlmostEqual(BASELINE,.0030676,places=7)

if __name__=="__main__":
    unittest.main()
