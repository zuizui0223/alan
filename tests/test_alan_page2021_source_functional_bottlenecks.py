"""Exact source-data algebra: pollinator q×V function is not visitation."""
import unittest
from scripts.analysis.alan_page2021_source_functional_bottlenecks_20261118 import (
 measures,aggregate,study_bootstrap,SCHEMA
)

class FunctionalBottleneckTests(unittest.TestCase):
    def test_function_narrower_when_common_visitors_have_no_pollen(self):
        # Common visitor V=90 but source q=0; rare V=10, q=50.
        r=measures([(90,0),(10,50)])
        self.assertAlmostEqual(r["functional_compositional_total_variation"],.9)
        self.assertAlmostEqual(r["visitation_effective_taxa_Hill2"],1/(.9**2+.1**2))
        self.assertAlmostEqual(r["pollen_delivery_effective_taxa_Hill2"],1.)
        self.assertAlmostEqual(r["top_visit_proxy_absolute_pollen_share_error"],.9)
        self.assertAlmostEqual(r["top_pollen_contributor_extra_share_over_top_visitor"],1.)
        self.assertTrue(r["functional_network_narrower"])

    def test_equal_efficiencies_do_not_generate_fake_functional_loss(self):
        a=measures([(70,4),(20,4),(10,4)])
        self.assertAlmostEqual(a["functional_compositional_total_variation"],0.,places=12)
        self.assertAlmostEqual(a["top_visit_proxy_absolute_pollen_share_error"],0.,places=12)
        self.assertAlmostEqual(a["top_pollen_contributor_extra_share_over_top_visitor"],0.,places=12)
        self.assertFalse(a["functional_network_narrower"])

    def test_source_zero_total_effectiveness_fails_closed(self):
        with self.assertRaises(ValueError):
            measures([(10,0),(20,0)])
        with self.assertRaises(ValueError):
            measures([(10,-2),(20,4)])
        with self.assertRaises(ValueError):
            measures([(10,1)])

    def test_descriptive_study_cluster_needs_five_independent_sources(self):
        r=[{"source_study":"A"}]*10
        self.assertEqual(study_bootstrap(r)["status"],
                         "fewer_than_five_source_studies_STOP")
        self.assertEqual(SCHEMA,
              "alan_page2021_33study_72plant_network_functional_tv_and_loss_v1")

    def test_group_summary_hill_and_functional_loss_gap(self):
        a=measures([(90,0),(10,50)])
        b={
            "source_study":"A","source_plant":"P","functional_total_variation":
                a["functional_compositional_total_variation"],
            "functional_Hill2":a["pollen_delivery_effective_taxa_Hill2"],
            "visitor_Hill2":a["visitation_effective_taxa_Hill2"],
            "functionally_narrower":a["functional_network_narrower"],
            "frequent_visitor_proxy_abs_error":
                a["top_visit_proxy_absolute_pollen_share_error"],
            "functional_top_taxon_pollen_share_gap":
                a["top_pollen_contributor_extra_share_over_top_visitor"],
            "most_visited_differs_most_functional":True,
            "original_rare_functional_taxon_present":False
        }
        s=aggregate([b])
        self.assertAlmostEqual(s["top_functional_share_gap_among_rank_mismatch_only"],1.)
        self.assertAlmostEqual(s["frequent_visitor_proxy_abs_pollen_share_error_mean"],.9)
        self.assertEqual(s["source_plant_network_cells"],1)
if __name__=="__main__":
    unittest.main()
