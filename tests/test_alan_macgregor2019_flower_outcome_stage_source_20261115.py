import unittest
from scripts.analysis.alan_macgregor2019_flower_outcome_stage_source_20261115 import compute

class TestMacgregorRealSourceStage(unittest.TestCase):
 def test_true_count_balance_and_separate_quality_denominator(self):
  z=compute()
  a=z["published_primary_counts"];q=z["genuine_staged_descriptive_proportions"]
  self.assertEqual(a["seed_present_flowers"]-a["larval_infested_seed_present_flowers"],a["uninfested_seed_capsules_quantified"])
  self.assertAlmostEqual(q["larval_infestation_per_seed_positive_flower"],46/194)
  self.assertAlmostEqual(q["seed_count_measured_per_total_mature_flower"],148/417)
  self.assertTrue(z["strong_larval_source_warning"]["seed_count_for_infested_capsules_MISSING_not_zero"])
 def test_published_diel_union_is_not_statistical_interaction(self):
  z=compute()["published_original_diel_pollination_success"]
  self.assertAlmostEqual(z["model_independent_union_from_day_and_night"],.705472)
  self.assertAlmostEqual(z["observed_open_minus_independence_union"],.014528)
  self.assertAlmostEqual(z["marginal_night_access_loss_1_minus_day_over_open"],1-.376/.720)
  self.assertTrue(z["independence_union_NOT_a_statistical_test"])
 def test_larval_fungibility_and_identity_unsupported(self):
  z=compute()
  p=z["published_primary_results"];l=z["strong_larval_source_warning"]
  self.assertFalse(p["larvae_confirmed_to_be_adult_Hadena_species"])
  self.assertEqual(p["one_raised_larva_identified_as"],"Lacanobia oleracea")
  self.assertTrue(l["not_a_measured_field_light_effect_on_seed_predation"])
  self.assertEqual(z["P1_x_P3_new_original_interaction_estimates"],0)
 def test_single_experiment_cannot_be_split_into_417_lamp_trials(self):
  z=compute()["replication_caveat"]
  self.assertTrue(z["flowers_not_independent_light_randomizations"])
  self.assertTrue(z["female_plants_reused_between_experimental_runs"])
if __name__=="__main__":unittest.main()
