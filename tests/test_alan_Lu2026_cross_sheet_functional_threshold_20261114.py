"""Mathematical tests for the genuine-source Lu cross-sheet audit."""
import unittest
from scripts.analysis.alan_Lu2026_cross_sheet_functional_threshold_20261114 import (
    analyze, cohort_sensitivity, weighted_two_order_ratio
)

def reported_aggregate_fixtures():
    def source(m0,m1,f0,f1,t0,t1):
        return {"Lepidoptera":{"original_dark_total_visit_events":m0,
                               "original_lit_total_visit_events":m1},
                "Diptera":{"original_dark_total_visit_events":f0,
                           "original_lit_total_visit_events":f1},
                "Total":{"original_dark_total_visit_events":t0,
                         "original_lit_total_visit_events":t1}}
    def plant(n,m0,m1,f0,f1,t0,t1):
        def item(a,b):
            return {"unlit_source_positive_visit_counts":a,
                    "lit_source_positive_visit_counts":b}
        return {"plant_taxa_with_both_lit_unlit_positive_visit_records":n,
                "treatment_visit_counts_UNADJUSTED_FOR_TRUE_FLOWER_AVAILABILITY":{
                    "moth":item(m0,m1),"fly":item(f0,f1),"raw_visits":item(t0,t1)}}
    survey={"Nocturnal":{"experimental_survey_pair_group_effects":
              source(838,357,54,151,1312,848)}}
    species={"by_period":{"Nocturnal":{
        "realized_control_moth_visitation_majority_groups":{
            "source_unlit_moth_fraction_above_half":plant(5,680,255,87,133,938,517),
            "source_unlit_moth_fraction_not_over_half":plant(6,95,73,19,48,249,231)},
        "cross_block_control_only_visitation_syndrome_transfer":{
            "n_supported_heldout_blocks":4,
            "heldout_blocks_where_moth_majority_stronger_relative_visitor_loss":2}}}}
    return survey,species

class LuOriginalCrossSheetAudit(unittest.TestCase):
    def test_raw_reported_aggregate_cross_sheet_incompatibility(self):
        survey,species=reported_aggregate_fixtures()
        out=analyze(survey,species)
        c=out["cross_sheet"]
        self.assertTrue(c["focal_11_plant_sum_exceeds_survey_total"])
        self.assertEqual(c["focal_11_plant_sums"]["fly"]["dark"],106)
        self.assertEqual(c["focal_11_plant_sums"]["fly"]["lit"],181)
        self.assertEqual(c["all_400_survey_counts"]["fly"]["dark"],54)
        self.assertEqual(c["all_400_survey_counts"]["fly"]["lit"],151)
        self.assertEqual(len(c["strict_same_unit_focal_subset_contradiction"]),2)
        self.assertFalse(out["generalization_gate"]["pooled_5_vs_6_classification_transfers_reliably"])

    def test_threshold_is_exact_no_pollen_claim(self):
        s=analyze(*reported_aggregate_fixtures())["within_sheet_visit_only_sensitivity"]
        expected={"400_survey_nocturnal":481/97,
                  "1038_plant_rows_focal_5_moth_majority":425/46,
                  "1038_plant_rows_focal_6_other":22/29}
        for cohort,value in expected.items():
            item=s[cohort]
            self.assertAlmostEqual(
              item["relative_fly_to_lepidoptera_per_visit_effectiveness_break_even_two_orders_only"],
              value)
            self.assertAlmostEqual(weighted_two_order_ratio(item["source_counts"],value),1.0)
            self.assertFalse(item["effective_pollen_or_fruit_measured"])

    def test_efficiency_scenario_only(self):
        s=analyze(*reported_aggregate_fixtures())["within_sheet_visit_only_sensitivity"]
        all_source=s["400_survey_nocturnal"]["source_counts"]
        self.assertAlmostEqual(weighted_two_order_ratio(all_source,1),508/892)
        self.assertLess(weighted_two_order_ratio(all_source,2),1)
        self.assertGreater(weighted_two_order_ratio(all_source,5),1)
        with self.assertRaises(ValueError):
            weighted_two_order_ratio(all_source,-1)

if __name__=="__main__":
    unittest.main()
