import importlib.util
import json
from pathlib import Path
import unittest

src=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_A_strict_readiness_20261115.py"
spec=importlib.util.spec_from_file_location("alanreadiness",src)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestALANReadiness(unittest.TestCase):
    def setUp(self):
        self.registry=json.loads(m.REG.read_text(encoding="utf-8"))
    def test_reproductive_publications_not_experiments(self):
        x=m.derive(self.registry)
        self.assertEqual(x["candidate_publications_total"],14)
        self.assertEqual(x["candidate_reproduction_papers"],12)
        self.assertEqual(x["confirmed_direct_lighting_reproduction_publications"],9)
        self.assertEqual(x["known_light_by_daynight_access_factorial_papers"],2)
    def test_lu_open_not_assumed_night(self):
        x=m.derive(self.registry)
        lu=next(r for r in x["source_statuses"] if r["source_id"]=="lu2026jpe")
        self.assertIn("NIGHT ONLY",lu["outcome_grain"])
        self.assertEqual(lu["qualified_primary_repro_lnRR_and_cluster_SE"],"no")
        self.assertIsNone(x["global_meta_beta_V"])
    def test_macgregor_not_dark_specific_v(self):
        x=m.derive(self.registry)
        mac=next(r for r in x["source_statuses"] if r["source_id"]=="macgregor2019")
        self.assertEqual(mac["dark_specific_V_and_SE_recovered"],"no")
    def test_changed_source_data_fails_closed(self):
        r=json.loads(json.dumps(self.registry))
        lu=next(x for x in r["entries"] if x["id"]=="lu2026jpe")
        lu["source_audit"]["Pedicularis_seed_lit_dark_open"] = .44
        with self.assertRaises(AssertionError):m.derive(r)
    def test_old_observational_not_promoted(self):
        r=json.loads(json.dumps(self.registry))
        bar=next(x for x in r["entries"] if x["id"]=="bariles2021")
        bar["light_design"]="manipulated"
        with self.assertRaises(AssertionError):m.derive(r)
    def test_cannot_relabel_mechanism_only_to_reproduction(self):
        r=json.loads(json.dumps(self.registry))
        ph=next(x for x in r["entries"] if x["id"]=="grognuz2026pheno")
        ph["a_primary_stage"]="A_endpoint_candidate"
        with self.assertRaises(AssertionError):m.derive(r)
    def test_confirmed_meta_beta_unavailable(self):
        x=m.derive(self.registry)
        self.assertEqual(x["open_reproduction_lnRR_with_clustered_SE_and_independent_V"],0)
        self.assertEqual(x["blockers_by_direct_study"]["dark_only_moderator_missing"],9)
if __name__=="__main__":unittest.main()
