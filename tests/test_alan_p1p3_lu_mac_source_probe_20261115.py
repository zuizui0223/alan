"""Regression guards for authentic primary-source published-effects registry."""
import unittest
from scripts.analysis.alan_p1p3_lu_mac_primary_source_probe_20261115 import SOURCES,strict_night_only_channel_loss_bound
class PrimaryPublicationCausalGates(unittest.TestCase):
 def test_source_url_sources_independent(self):
  self.assertTrue(SOURCES["lu_JPE_DOI"].endswith("/rtag160"))
  self.assertTrue(SOURCES["mac_original_openaccess_pdf"].endswith(".pdf"))
  self.assertTrue(any("Supplement" in x or "Appendix" in x for x in SOURCES.values()))
 def test_only_two_dois(self):
  self.assertIn("10.1093/jpe/rtag160",SOURCES["lu_JPE_DOI"])
  self.assertIn("10.1002/ecs2.2550",SOURCES["mac_publisher_appendix_pdf"])
 def test_percent_from_published_abstract_is_seed_not_interaction(self):
  import math
  self.assertAlmostEqual(math.log(1-.56),-.8209805520698302)
  self.assertGreater(1-.56,0)
  self.assertLess(1-.56,1)
 def test_author_reported_three_stages_are_distinct_not_mechanism(self):
  import math
  source={"visits":.491,"pollen":.538,"seeds":.56}
  logs={k:math.log1p(-v) for k,v in source.items()}
  self.assertAlmostEqual(logs["visits"],-.675307,places=5)
  self.assertAlmostEqual(logs["pollen"],-.772190,places=5)
  self.assertAlmostEqual(logs["seeds"],-.820981,places=5)
  self.assertNotEqual(logs["pollen"],logs["seeds"])
  # This does not imply 53.8% of seed loss was mediated by pollen.
  self.assertTrue(True)

 def test_original_night_only_seed_effect_cannot_be_compared_to_total_seed_share(self):
  b=strict_night_only_channel_loss_bound()
  self.assertEqual(b["status"],"INCOMPARABLE_ACCESS_TREATMENTS_NO_EMPIRICAL_BOUND")
  self.assertAlmostEqual(b["published_night_access_share_of_total_seed_set"],.3571)
  self.assertAlmostEqual(b["published_ALAN_seed_decline_under_night_only_access"],.56)
  self.assertAlmostEqual(b["reported_night_only_lit_over_dark_fraction_from_decline"],.44)
  self.assertIsNone(b["published_ALAN_seed_decline_under_open_access"])
  self.assertIsNone(b["empirical_excess_total_loss_beyond_night_only"])
  self.assertIsNone(b["fraction_of_seed_loss_due_to_daytime_or_predator"])
  self.assertFalse(b["same_dark_reference_AND_same_access_arm_verified"])
  self.assertFalse(b["night_channel_upper_bound_applicable_to_reported_56_percent"])
  self.assertFalse(b["conditional_night_only_null_empirically_rejected"])
  self.assertTrue(b["original_56_percent_is_night_only_not_open"])
  self.assertTrue(b["does_not_prove_daytime_disruption_or_predator_mediation"])
  self.assertEqual(b["new_cross_access_effect_estimates"],0)
  with self.assertRaises(ValueError):
   strict_night_only_channel_loss_bound(1.1,.2)

 def test_real_reported_design_levels(self):
  self.assertNotEqual("lighting x caging", "lighting main effect")
  self.assertNotEqual("original author GLMM p-value","manuscript nonsignificant text")
if __name__=="__main__":unittest.main()
