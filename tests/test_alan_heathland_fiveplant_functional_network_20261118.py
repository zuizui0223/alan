"""Fail-closed original five-plant functional pollen concentration mechanics."""
import unittest
from scripts.analysis.alan_heathland_fiveplant_functional_pollen_network_20261118 import (
    summarize_plant,day_cluster_bootstrap,SCHEMA
)
class HeathNetworkUnitTests(unittest.TestCase):
    def make(self,plant,visitors):
        return [{"plant":plant,"taxon":taxon,"guild":guild,
                 "date":day,"conspecific_pollen":float(pollen)}
                for taxon,guild,day,pollen in visitors]
    def test_fewer_visits_can_contribute_more_pollen(self):
        rows=self.make("Erica_tetralix",[
          ("Apis_mellifera","bees","d1",1),
          ("Apis_mellifera","bees","d2",1),
          ("Apis_mellifera","bees","d3",1),
          ("Bombus_jonellus","bees","d2",10),
        ])
        r=summarize_plant(rows,"Erica_tetralix")
        self.assertFalse(r["source_same_taxon_tops_both"])
        self.assertEqual(r["top_source_visitation_taxon"],"Apis_mellifera")
        self.assertEqual(r["top_source_functional_pollen_depositor"],"Bombus_jonellus")
        self.assertAlmostEqual(r["top_visitor_functional_pollen_fraction"],3/13)
        self.assertAlmostEqual(r["top_functional_taxon_pollen_fraction"],10/13)

    def test_rare_original_visitor_detected_by_functional_share(self):
        rows=self.make("Erica_tetralix",
            [("Apis_mellifera","bees",f"d{i}",1) for i in range(20)]+
            [("Bombus_jonellus","bees","d0",10)])
        r=summarize_plant(rows,"Erica_tetralix")
        rare=r["original_rare_under5pct_visits_over10pct_pollen_taxa"]
        self.assertEqual(len(rare),1)
        self.assertEqual(rare[0]["source_taxon"],"Bombus_jonellus")
        self.assertAlmostEqual(rare[0]["source_visit_fraction"],1/21)
        self.assertAlmostEqual(rare[0]["source_pollen_fraction"],1/3)

    def test_with_no_fly_records_no_fake_zero_efficiency(self):
        rows=self.make("Erica_tetralix",[
           ("Apis_mellifera","bees",f"d{i}",2) for i in range(10)])
        r=summarize_plant(rows,"Erica_tetralix")
        self.assertEqual(r["hypothetical_10pct_bee_visits_replaced_with_fly_visits"]["status"],
                         "source_fly_within_plant_replication_STOP")

    def test_functional_fly_swap_is_explicitly_counterfactual(self):
        rows=self.make("Calluna_vulgaris",
          [("Apis_mellifera","bees",f"d{i}",10) for i in range(12)]+
          [("muscid","flies",f"d{i%6}",5) for i in range(12)])
        r=summarize_plant(rows,"Calluna_vulgaris")
        q=r["hypothetical_10pct_bee_visits_replaced_with_fly_visits"]
        self.assertEqual(q["status"],"source_support_pass_conditional_only")
        self.assertAlmostEqual(q["hypothetical_change_in_bee_fly_contributed_pollen_fraction"],-.05)
        self.assertTrue(q["not_moth_or_ALAN_pollinator_replacement"])

if __name__=="__main__":
    unittest.main()
