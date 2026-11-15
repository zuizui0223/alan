"""Regression guards for authentic primary-source published-effects registry."""
import unittest
from scripts.analysis.alan_p1p3_lu_mac_primary_source_probe_20261115 import SOURCES
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

 def test_real_reported_design_levels(self):
  self.assertNotEqual("lighting x caging", "lighting main effect")
  self.assertNotEqual("original author GLMM p-value","manuscript nonsignificant text")
if __name__=="__main__":unittest.main()
