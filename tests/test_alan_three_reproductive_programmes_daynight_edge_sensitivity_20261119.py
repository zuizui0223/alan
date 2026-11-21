import copy
import json
import math
import unittest
from pathlib import Path
from scripts.analysis.alan_three_reproductive_programmes_daynight_edge_sensitivity_20261119 import report

LEDGER=Path("results/all_species_forest/biotic_interactions/alan_three_true_plant_reproduction_published_coefficient_ledger_20261119.json")


class OriginalSileneAndYuccaReproductionPathwayTests(unittest.TestCase):
    def setUp(self):
        self.data=json.loads(LEDGER.read_text())

    def test_original_published_caging_complementarity_is_only_scenario(self):
        x=report(self.data)
        self.assertEqual(x["primary_published_programmes_n"],3)
        m=x["macgregor2019_published_day_night_functional_pathways"]
        expected=1-(1-.376)*(1-.528)
        self.assertAlmostEqual(m["independent_day_night_route_predicted_open"],expected)
        self.assertAlmostEqual(m["observed_minus_independent_route_prediction_percentage_points"],100*(.72-expected))
        self.assertAlmostEqual(m["published_damaged_fraction_among_successfully_pollinated"],46/194)
        self.assertTrue(m["marginal_caging_rates_not_a_dark_control_only_comparison"])
        s=m["night_loss_sensitivity_not_a_causal_ALAN_dose_response"]
        self.assertEqual(len(s),5)
        self.assertAlmostEqual(s[2]["predicted_combined_pollination_success_if_day_unchanged"],1-(1-.376)*(1-.264))
        self.assertAlmostEqual(s[-1]["predicted_combined_pollination_success_if_day_unchanged"],.376)
        self.assertTrue(all(row["not_real_ALAN_treatment_effect"] for row in s))

    def test_Giavi_published_model_effect_stage_not_illegal_seed_conversion(self):
        d=report(self.data)
        g=d["giavi2020_original_field_and_lab_model_table_contrasts"]
        self.assertEqual(g["physical_lit_meadows"],4)
        self.assertEqual(g["physical_dark_control_meadows"],4)
        q=g["original_published_effects_separately_by_endpoint"]
        self.assertEqual(len(q),7)
        v=next(row for row in q if row["endpoint"]=="Silene_infested_fruit_fraction" and row["orientation"]=="adjacent/dark")
        self.assertAlmostEqual(v["source_model_conditional_attack_odds_ratio"],math.exp(1.90))
        self.assertAlmostEqual(v["attack_odds_ratio_Wald_95pct_CI"][0],math.exp(1.90-1.96*.576))
        self.assertAlmostEqual(v["attack_odds_ratio_Wald_95pct_CI"][1],math.exp(1.90+1.96*.576))
        x=next(row for row in q if row["endpoint"]=="Silene_uninfected_capsules_per_exposed_flower" and row["orientation"]=="dark/adjacent")
        self.assertAlmostEqual(x["published_beta"],1.04)
        self.assertEqual(x["published_p_value"],.037)
        self.assertTrue(x["do_not_exponentiate_as_fruit_or_seed_probability_ratio"])
        self.assertNotIn("source_model_conditional_attack_odds_ratio",x)
        self.assertTrue(g["cannot_pretend_published_reproduction_link_is_actual_seed_loss_ratio"])

        # Model-based necessary condition; dark-site infection rates here are
        # hypothetical and cannot be confused with genuine Giavi plot rates.
        thresholds=g["conditional_compensatory_pollination_thresholds"]
        self.assertEqual(len(thresholds),4)
        t=next(row for row in thresholds
               if row["hypothetical_original_dark_fruit_infection_probability"]==.10)
        OR=math.exp(1.9)
        self.assertAlmostEqual(
            t["minimum_adjacent_vs_dark_capsule_formation_probability_ratio_for_no_uninfected_capsule_loss"],
            1+.1*(OR-1))
        self.assertAlmostEqual(
            t["model_constrained_adjacent_fruit_infection_probability"],
            OR*.1/(.9+OR*.1))
        self.assertAlmostEqual(
            t["maximum_original_dark_capsule_formation_probability_where_full_rescue_is_even_possible"],
            1/(1+.1*(OR-1)))
        self.assertTrue(all(row["NOT_an_observed_ALAN_pollination_effect"] for row in thresholds))
        self.assertTrue(all(row["baseline_infection_hypothetical_not_from_Giavi_original_samples"] for row in thresholds))

        self.assertFalse(d["primary_reproductive_effect_random_effects_meta_fitted"])

    def test_Wilson_skyglow_not_randomized_lamp_or_original_effect_SE(self):
        r=report(self.data)
        w=r["wilson2021_distinguish_ambient_skyglow_and_direct_tower"]
        self.assertEqual(len(w),4)
        self.assertEqual(len([x for x in w if x["exposure"]=="skyglow"]),2)
        self.assertTrue(all(x["source_light_assignment"] is False for x in w if x["exposure"]=="skyglow"))
        self.assertTrue(all(x["source_effect_size_beta_SE_not_verified"] for x in w))
        self.assertEqual(r["three_primary_effects_with_same_endpoint_SE_and_verified_treatment_units"],0)

    def test_never_derive_fake_treatment_counts_or_SE_from_caging(self):
        data=copy.deepcopy(self.data)
        m=next(x for x in data["primary_studies"] if x["programme"]=="Macgregor2019")
        m["design"]["source_original_numerical_treatment_CSV_available"]=True
        with self.assertRaisesRegex(ValueError,"missing original treatment"):
            report(data)
        data=copy.deepcopy(self.data)
        g=next(x for x in data["primary_studies"] if x["programme"]=="Giavi2020")
        g["design"]["adjacent_dark_shares_treatment_meadow"]=False
        with self.assertRaisesRegex(ValueError,"separate randomized"):
            report(data)
        data=copy.deepcopy(self.data)
        g=next(x for x in data["primary_studies"] if x["programme"]=="Giavi2020")
        g["published_model_contrasts"][1]["SE"]=0.
        with self.assertRaisesRegex(ValueError,"beta or SE"):
            report(data)

if __name__=="__main__":
    unittest.main()
