import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
REG=json.loads((ROOT/"config/alan_A_reproductive_screen_20261115.json").read_text())
ENTRIES=REG["entries"]

class TestWilsonLightGrain(unittest.TestCase):
    def test_registry_counts(self):
        self.assertEqual(REG["reviewed_reproductive_papers"],12)
        self.assertEqual(REG["design_counts_among_reproductive"]["manipulated"],9)
        self.assertEqual(REG["design_counts_among_reproductive"]["observational"],3)
    def test_one_wilson_publication(self):
        pubs=[r for r in ENTRIES if r["doi"]=="10.1093/icb/icab010"]
        self.assertEqual(len(pubs),1)
        self.assertTrue("manipulated" in pubs[0]["light_design"])
    def test_true_independent_repro_layout(self):
        wilson=next(r for r in ENTRIES if r["id"]=="wilson2021")
        a=wilson["source_audit"]
        self.assertEqual(a["n_manipulated_transects"],a["n_lit_transects"]+a["n_dark_transects"])
        self.assertEqual(a["n_manipulated_transects"],15)
        self.assertEqual(a["n_field_sites"],5)
        self.assertEqual(a["n_fruit_response_plants"],58)
        self.assertFalse(a["skyglow_randomized"])
        self.assertTrue(a["exclude_light_trap_data_from_repro_effect"])
    def test_no_mislabel_85pct_and_no_fake_moderator(self):
        wilson=next(r for r in ENTRIES if r["id"]=="wilson2021")
        self.assertEqual(wilson["source_audit"]["published_illuminance_CI85"],[.019,.657])
        self.assertFalse(wilson["verified_joint_effect_moderator_se"])
        self.assertEqual(wilson["external_night_dependence"],"none_proven")

if __name__=="__main__":
    unittest.main()
