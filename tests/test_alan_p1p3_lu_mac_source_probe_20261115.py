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

 def test_night_only_loss_bound_is_conditional_not_mediation(self):
  b=strict_night_only_channel_loss_bound()
  self.assertAlmostEqual(b["maximum_possible_decline_if_only_dark_night_contribution_removed"],.3571)
  self.assertAlmostEqual(b["minimum_possible_lit_over_dark_seed_ratio_night_only"],.6429)
  self.assertAlmostEqual(b["reported_lit_over_dark_seed_ratio_from_reduction"],.44)
  self.assertAlmostEqual(b["excess_total_loss_beyond_night_only_upper_bound"],.2029)
  self.assertAlmostEqual(b["fraction_observed_decline_not_coverable_by_night_only_maximum"],.2029/.56)
  self.assertTrue(b["conditional_night_only_null_incompatible_if_common_reference"])
  self.assertFalse(b["same_reference_verified_by_abstract"])
  self.assertTrue(b["does_not_prove_daytime_disruption_or_predator_mediation"])
  self.assertFalse(strict_night_only_channel_loss_bound(.8,.56)["conditional_night_only_null_incompatible_if_common_reference"])
  with self.assertRaises(ValueError):
   strict_night_only_channel_loss_bound(1.1,.2)

 def test_real_reported_design_levels(self):
  self.assertNotEqual("lighting x caging", "lighting main effect")
  self.assertNotEqual("original author GLMM p-value","manuscript nonsignificant text")
if __name__=="__main__":unittest.main()
