import copy
import json
import math
import unittest
from pathlib import Path

from scripts.analysis.alan_stageA_wetland_land_process_source_validity_20261119 import audit
SOURCE=Path("config/alan_stageA_water_land_process_specific_primary_gate_20261119.json")

class StageAAquaticTerrestrialSourceValidation(unittest.TestCase):
    def setUp(self):
        self.a=json.loads(SOURCE.read_text())
    def test_six_source_programmes_not_seven_independent_papers(self):
        d=audit(self.a)
        self.assertEqual(d["targeted_original_programme_k"],6)
        self.assertEqual(d["targeted_publication_n"],7)
        self.assertEqual(d["original_programme_habitat_k_not_a_meta_denominator"]["freshwater_riparian"],2)
        self.assertEqual(d["original_programme_habitat_k_not_a_meta_denominator"]["freshwater_headwater_stream"],1)
        self.assertEqual(d["original_programme_habitat_k_not_a_meta_denominator"]["marine_shallow_water"],1)
        self.assertEqual(d["original_programme_habitat_k_not_a_meta_denominator"]["terrestrial_forest"],2)
        self.assertFalse(d["independent_A_wetland_by_terrestrial_meta_moderator_fit"])
        self.assertTrue(d["no_meta_pooled_mean_CI_or_pvalue"])
    def test_lit_emergence_vs_flying_capture_are_distinct_not_variance_ready(self):
        d=audit(self.a)
        self.assertAlmostEqual(d["Westhavelland_2014_Diptera_emergence_CPUE_lit_dark_descriptive_mean_ratio"],3.72/7.84)
        self.assertAlmostEqual(d["Westhavelland_2014_Diptera_emergence_descriptive_log_mean_ratio_NO_INDEPENDENT_SITE_SE"],math.log(3.72/7.84))
        self.assertEqual(d["Westhavelland_2014_flying_Diptera_catch_approximately_light_dark_ratio"],8)
        self.assertNotEqual(d["Westhavelland_2014_Diptera_emergence_CPUE_lit_dark_descriptive_mean_ratio"],
                            d["Westhavelland_2014_flying_Diptera_catch_approximately_light_dark_ratio"])
        self.assertTrue(d["no_sample_size_from_nested_streetlamps_or_traps"])
    def test_source_site_pair_and_four_streams_not_lamp_replication(self):
        x=copy.deepcopy(self.a)
        w=next(z for z in x["programmes"] if z["programme_id"]=="WesthavellandDrainageDitch2012")
        w["original_n_treated_sites"]=12
        with self.assertRaisesRegex(ValueError,"one lit and one dark"):
            audit(x)
        x=copy.deepcopy(self.a)
        s=next(z for z in x["programmes"] if z["programme_id"]=="Perkin2014_4paired_headwater_streams")
        s["original_n_independent_streams"]=8
        with self.assertRaisesRegex(ValueError,"four paired streams"):
            audit(x)
        x=copy.deepcopy(self.a)
        m=next(z for z in x["programmes"] if z["programme_id"].startswith("NavarroBarranco"))
        m["site_type"]="freshwater_riparian"
        with self.assertRaisesRegex(ValueError,"marine amphipod"):
            audit(x)
    def test_predator_count_not_plant_seed_mediation(self):
        x=audit(self.a)
        self.assertTrue(x["no_habitat_sign_claim_without_matched_endpoint_and_independent_site_replicates"])
        self.assertTrue(x["marine_amphipod_lighttrap_attraction_not_freshwater_emergence"])
        self.assertEqual(x["Perkin2014_emergence_and_trout_growth"],"not resolved in published source, NOT zero")

if __name__=="__main__":
    unittest.main()
