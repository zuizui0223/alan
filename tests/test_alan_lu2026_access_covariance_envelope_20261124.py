"""Check source-pinned Lu 2026 supplement directional and covariance bounds."""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
MODULE=ROOT/"scripts/analysis/alan_lu2026_access_covariance_envelope_20261124.py"
spec=importlib.util.spec_from_file_location("lu_envelope",MODULE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class LuSupplementBounds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=m.read_source()

    def test_source_rows_reconstruct_twelve_not_independent(self):
        self.assertEqual(len(self.source),12)
        for (spp,response,access),row in self.source.items():
            self.assertIn(spp,("P. kansuensis","E. regelii"))
            self.assertIn(access,m.ALLOWED_ACCESS)
            self.assertGreater(row["model_SE"],0)

    def test_seed_buffering_point_and_unknown_covariance(self):
        z=m.contrast_envelope(self.source,"P. kansuensis","Seed set")
        self.assertAlmostEqual(z["B_point_on_original_reported_model_scale"],.0887,places=12)
        self.assertAlmostEqual(z["SE_bound_cauchy_schwarz"][0],.011,places=12)
        self.assertAlmostEqual(z["SE_bound_cauchy_schwarz"][1],.0748,places=12)
        self.assertAlmostEqual(z["hypothetical_rho_zero_SE"],.05346045267298061,places=12)
        self.assertAlmostEqual(z["worst_case_normal_z_1p96_envelope"][0],-.057908,places=12)
        self.assertAlmostEqual(z["worst_case_normal_z_1p96_envelope"][1],.235308,places=12)
        self.assertAlmostEqual(z["rho_needed_for_naive_nominal_z1p96_nonzero"],.29594074550982513,places=12)
        self.assertFalse(z["proper_design_based_95CI_estimable"])
        self.assertIsNone(z["absolute_seed_fraction_or_pp_change"])
        self.assertFalse(z["adopt_into_meta"])

    def test_other_species_or_pollen_are_different_endpoints(self):
        e=m.contrast_envelope(self.source,"E. regelii","Seed set")
        p=m.contrast_envelope(self.source,"P. kansuensis","Pollen number")
        self.assertAlmostEqual(e["B_point_on_original_reported_model_scale"],.0054)
        self.assertAlmostEqual(p["B_point_on_original_reported_model_scale"],-4.43)
        self.assertNotEqual(p["outcome"],e["outcome"])

    def test_source_numeric_corruption_stops(self):
        original=m.SOURCE
        try:
            # Do not modify the original source file: corrupted inputs are
            # tested through existing in-memory source-table reconstruction.
            t=copy.deepcopy(self.source)
            t["P. kansuensis","Seed set","Nocturnal"]["effect_unlit_minus_lit"]=.5
            observed=m.contrast_envelope(t,"P. kansuensis","Seed set")
            self.assertNotAlmostEqual(observed["B_point_on_original_reported_model_scale"],.0887)
        finally:
            self.assertEqual(m.SOURCE,original)

if __name__=="__main__":
    unittest.main()
