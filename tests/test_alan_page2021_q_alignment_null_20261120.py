"""Synthetic q-V conditional-null invariants; no field-ALAN claims."""
import random
import unittest
from scripts.analysis.alan_page2021_visitor_q_alignment_null_20261120 import (
 thresholds, permute, network_null, summarize, study_bootstrap
)
class QVisitorNull(unittest.TestCase):
    def test_unique_highest_q_on_most_frequent_is_best_failure(self):
        self.assertEqual(thresholds([12.,2.,3.],[8.,1.,2.]),(1,1))
        self.assertEqual(thresholds([12.,2.,3.],[1.,8.,2.]),(0,0))

    def test_efficiency_ties_cannot_trigger_best_failure(self):
        b,p=thresholds([10.,2.,1.],[4.,4.,1.])
        self.assertEqual(b,0)

    def test_permutation_preserves_empirical_efficiency_multiset(self):
        q=[10.,1.,3.,4.]
        g=["fly","fly","bee","moth"]
        for restricted in (False,True):
            rng=random.Random(13)
            out=permute(q,g,rng,restricted)
            self.assertEqual(sorted(q),sorted(out))
            if restricted:
                self.assertEqual(q[2:],out[2:])

    def test_unique_max_exact_probability_one_over_taxa(self):
        data={"visits":[10.,6.,2.],"q":[9.,3.,1.],
              "guilds":["bee","fly","moth"],
              "source_study":"A","source_plant":"Plant"}
        o=network_null(data,random.Random(1),reps=5997,restricted=False)
        self.assertAlmostEqual(o["conditional_best_failure_null_exact_if_unique_max_q"],1/3)
        self.assertLess(abs(o["null_probability_best_failure"]-1/3),.03)
        self.assertEqual(o["observed_best_failure"],1)

    def test_guild_restricted_no_change_when_all_singletons(self):
        data={"visits":[10.,6.,2.],"q":[9.,3.,1.],
              "guilds":["bee","fly","moth"],
              "source_study":"A","source_plant":"Plant"}
        o=network_null(data,random.Random(1),reps=399,restricted=True)
        self.assertEqual(o["null_replicates_changed_q_assignment"],0)
        self.assertEqual(o["null_probability_best_failure"],1.0)

    def test_study_cluster_is_not_taxon_repeat(self):
        data=[{"source_study":"S1","observed_best_failure":1,
               "observed_proportional_failure":1,
               "null_probability_best_failure":.3,
               "null_probability_proportional_failure":.5},
              {"source_study":"S1","observed_best_failure":1,
               "observed_proportional_failure":1,
               "null_probability_best_failure":.4,
               "null_probability_proportional_failure":.6},
              {"source_study":"S2","observed_best_failure":0,
               "observed_proportional_failure":0,
               "null_probability_best_failure":.2,
               "null_probability_proportional_failure":.4}]
        agg,studies=summarize(data)
        self.assertEqual(agg["networks"],3)
        self.assertEqual(agg["source_studies"],2)
        self.assertEqual(len(studies),2)
        intervals=study_bootstrap(studies,299)
        self.assertEqual(len(intervals["best_observed_minus_null"]),2)

if __name__=="__main__":
    unittest.main()
