"""Source-locked wet flora food web coverage: absence and pseudoreplication must not be invented."""
import unittest
from pathlib import Path
from scripts.analysis.alan_wet_biotic_multiguild_source_20261118 import (
    FILES,EXPECTED_GUILDS,audit
)
class TestWetFoodwebDocumentaryGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a=audit(Path("results/all_species_forest"))

    def test_frozen_wet_plant_denominator(self):
        a=self.a["two_step_biotic_graph_support"]
        self.assertEqual(a["source_wet_flora_species"],363)
        self.assertEqual(a["source_candidate_plants"],58)
        self.assertEqual(a["source_other_wet_plants"],305)
        self.assertEqual(a["documented_source_matched_adult_visitor_predator_pairs"],272)
        self.assertEqual(sum(x["shared_visitor_predator_pair_types"] for x in self.a["predator_guilds"].values()),272)

    def test_non_inat_provenance_is_not_independent_study(self):
        x=self.a["two_step_biotic_graph_support"]
        self.assertEqual(x["chains_with_at_least_one_non_iNaturalist_source_pair"],2)
        self.assertEqual(x["chains_with_only_iNaturalist_containing_source_pairs"],270)
        self.assertTrue(x["all_source_pairs_are_evidence_of_global_links_not_same_field_plot"])

    def test_adult_moth_spider_riparian_edges_do_not_establish_causal_chain(self):
        a=self.a["documented_riparian_moth_adult_spider_chains"]
        self.assertEqual(a["original_moth_to_spider_source_rows"],8)
        self.assertEqual(a["source_diet_edges_reliant_on_iNaturalist"],8)
        self.assertFalse(a["original_local_plant_spider_pollen_transfer_same_place_year_verified"])
        self.assertEqual(self.a["predator_guilds"]["spiders"]["shared_visitor_predator_pair_types"],91)
        self.assertEqual(self.a["predator_guilds"]["reptiles"]["shared_visitor_predator_pair_types"],3)

    def test_source_files_sha_locked(self):
        self.assertEqual(len(FILES),3)
        self.assertEqual(dict((g, x["shared_visitor_predator_pair_types"])
                              for g,x in self.a["predator_guilds"].items()),EXPECTED_GUILDS)

if __name__=="__main__":
    unittest.main()
