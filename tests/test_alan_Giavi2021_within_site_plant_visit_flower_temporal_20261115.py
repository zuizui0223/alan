"""Source synthetic fixtures; true original data only used in CI source workflow."""
import unittest
import pandas as pd
from scripts.analysis.alan_Giavi2021_within_site_plant_visit_flower_temporal_20261115 import (
    normalized_event_overlap, pair_plants, summarize_paired, site_bootstrap
)
class TemporalFunctionalTests(unittest.TestCase):
    def test_same_temporal_profiles_overlap_one_and_shift_zero(self):
        self.assertAlmostEqual(normalized_event_overlap([1,3,2],[2,6,4]),1.0)
        self.assertAlmostEqual(normalized_event_overlap([2,0,0],[0,0,5]),0.)
    def test_no_visits_means_unknown_not_zero_alignment(self):
        self.assertIsNone(normalized_event_overlap([2,3,1],[0,0,0]))
        self.assertIsNone(normalized_event_overlap([0,0,0],[1,2,3]))
        with self.assertRaises(ValueError):
            normalized_event_overlap([1,-2,3],[1,0,2])
    def test_common_plant_site_occasion_filter_prevents_imputed_missing(self):
        q=pd.DataFrame([
          {"site":"A","plant":"Test species","pair":i,
           "floral_cover_C":f,"floral_cover_T":f,
           "all3_visits_C":v,"all3_visits_T":t}
          for i,(f,v,t) in enumerate([(1,2,0),(4,4,1),(1,0,4)],1)])
        r=pair_plants(q,"all3_visits",3)
        self.assertEqual(len(r),1)
        self.assertEqual(r[0]["source_common_occasion_n"],3)
        self.assertGreater(r[0]["lit_overlap"],0)
        self.assertLess(r[0]["lit_overlap"],1)
        self.assertEqual(pair_plants(q,"all3_visits",4),[])
        self.assertAlmostEqual(normalized_event_overlap([1,4,1],[2,4,0]),.75)
    def test_only_site_resampling_and_small_n_gate(self):
        self.assertIsNone(site_bootstrap([1,2]))
        a=site_bootstrap([.1,.1,.1,.1,.1,.1],iterations=100)
        self.assertEqual(len(a),2)
        self.assertAlmostEqual(a[0],.1)
        x=summarize_paired([
          {"source_site":str(i),"plant":"P",
           "source_common_occasion_n":3,
           "lit_overlap":.5,
           "dark_overlap":.3,
           "lit_minus_dark_temporal_overlap":.2,
           "lit_total_visits_same_support":4,
           "dark_total_visits_same_support":3}
          for i in range(6)],11)
        self.assertEqual(x["supported_independent_sites"],6)
        self.assertTrue(x["not_a_pollinator_efficacy_or_true_phenological_mismatch"])
        self.assertTrue(x["source_missing_excluded_not_zero"])
if __name__=="__main__":unittest.main()
