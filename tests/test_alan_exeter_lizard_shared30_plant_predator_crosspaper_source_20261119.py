import copy
import json
import unittest
from pathlib import Path
import urllib.error

from scripts.analysis.alan_exeter_lizard_shared30_plant_predator_crosspaper_source_20261119 import (
    validate, inspect_html_links, source_scan
)

REG = Path("results/all_species_forest/biotic_interactions/alan_exeter_lizard_shared_30lightplots_crosspaper_provenance_20261119.json")


class SharedThirtyLightPlotsTests(unittest.TestCase):
    def setUp(self):
        self.registry=json.loads(REG.read_text())

    def test_identical_original_physical_light_sites_do_not_count_twice(self):
        result=validate(self.registry)
        self.assertEqual(result["field_programme"],"ExeterLizardLED30_2012_2015")
        self.assertEqual(result["published_article_count"],2)
        self.assertEqual(result["independently_randomized_light_field_experiments"],1)
        self.assertEqual(result["independently_randomized_physical_field_plots"],30)
        self.assertEqual(result["physical_plots_per_lighting_regime"],6)
        self.assertEqual(result["physical_year_overlap"],[2012,2013,2014])
        self.assertEqual(result["original_plant_species_with_flowering_dates_n"],3)
        self.assertEqual(result["original_spider_beetle_caught_individuals"],5180)
        self.assertEqual(result["original_spider_beetle_taxa"],136)
        self.assertFalse(result["original_plot_ID_codelist_usable_for_joint_source_model"])
        self.assertFalse(result["actual_source_joint_cross_trophic_paired_model_fitted"])
        self.assertTrue(result["no_predator_mediation_or_species_seed_effect_estimated"])

    def test_original_light_regime_aliases_preserved_between_two_papers(self):
        d=validate(self.registry)
        arms=d["light_arms"]
        self.assertEqual(set(arms),{"CONTROL","HIW","LIW","LIWT","AMB"})
        self.assertEqual(arms["LIW"]["Bennie2018"],"LIW")
        self.assertEqual(arms["LIW"]["Davies2017"],"DW")
        self.assertEqual(arms["LIWT"]["Bennie2018"],"LIWT")
        self.assertEqual(arms["LIWT"]["Davies2017"],"DWT")
        self.assertEqual(sum(r["physical_plots"] for r in arms.values()),30)

    def test_source_mutations_and_false_joint_fits_are_rejected(self):
        r=copy.deepcopy(self.registry)
        r["site"]["original_independently_randomized_light_plots"]=54
        with self.assertRaisesRegex(ValueError,"30"):
            validate(r)
        r=copy.deepcopy(self.registry)
        r["physical_treatment_arms"][2]["Davies2017"]="LIW"
        with self.assertRaisesRegex(ValueError,"aliases"):
            validate(r)
        r=copy.deepcopy(self.registry)
        r["source_linkage_evidence"]["independent_light_experiment_count"]=2
        with self.assertRaisesRegex(ValueError,"one original"):
            validate(r)
        r=copy.deepcopy(self.registry)
        r["new_joint_cross_trophic_estimate_gate"]["independent_crosspaper_same_plot_match_performed"]=True
        with self.assertRaisesRegex(ValueError,"before claiming"):
            validate(r)

    def test_archive_html_link_parser_accepts_official_exeter_only(self):
        html=b"""
        <a href='/repository/bitstream/handle/10871/27138/legitimate_original.csv'>CSV</a>
        <a href='https://ore.exeter.ac.uk/repository/handle/10871/25044'>Record</a>
        <a href='https://example.org/bitstream/handle/10871/27138/fake.csv'>Fake</a>
        <a href='/repository/bitstream/handle/10871/25044/OriginalData.RData'>RData</a>
        <a href='/repository/bitstream/handle/10871/27138/legitimate_original.csv'>Duplicate</a>
        """
        candidates=inspect_html_links(html,"https://ore.exeter.ac.uk/repository/handle/10871/27138")
        self.assertEqual(len(candidates),2)
        self.assertTrue(all(x.startswith("https://ore.exeter.ac.uk/") for x in candidates))
        self.assertTrue(any(x.endswith("legitimate_original.csv") for x in candidates))
        self.assertTrue(any(x.endswith("OriginalData.RData") for x in candidates))

    def test_source_access_status_does_not_claim_plot_level_join(self):
        def reader(url):
            if url.endswith("/10871/27138"):
                return b"<a href='/repository/bitstream/10871/27138/original.xlsx'>Dataset</a>",url,"text/html"
            raise urllib.error.HTTPError(url,404,"not found",None,None)
        d=source_scan(reader=reader)
        self.assertEqual(d["physical_independently_assigned_light_plots"],30)
        self.assertFalse(d["joint_original_plant_arachnid_plot_records_verified"])
        self.assertFalse(d["verified_source_joint_multi_trophic_treatment_effect_estimated"])
        self.assertTrue(d["no_original_plot_ID_or_arthropod_row_exported"])
        success=d["public_archival_page_access"]["Bennie2018_flowering_original_archived_handle"][0]
        self.assertEqual(success["status"],"SOURCE_METADATA_RETRIEVED")
        self.assertEqual(len(success["official_archived_file_candidates"]),1)
        self.assertTrue(d["metadata_success_not_raw_multitrophic_join"])


if __name__=="__main__":
    unittest.main()
