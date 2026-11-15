import importlib.util
import json
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_Silene_lat_dark_only_V_gate_20261115.py"
spec=importlib.util.spec_from_file_location("silenedark",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestSileneV(unittest.TestCase):
    def setUp(self):
        self.source=m.load_rows(m.ORIG)
        self.original=m.load_rows(m.SPEC)
    def test_original_42_source_rows(self):
        self.assertEqual(len(self.source),42)
        self.assertEqual(len(self.original),12)
    def test_baseline_all_studies_reproduced(self):
        r=m.compute(self.source,self.original)
        self.assertEqual(r["species_n"],12)
        self.assertAlmostEqual(r["s_lat_original_V"],.6718269901631905)
    def test_specific_lighting_contamination_repaired(self):
        r=m.compute(self.source,self.original)
        self.assertAlmostEqual(r["s_lat_strict_V"],.7236903716106219)
        self.assertAlmostEqual(r["s_lat_absolute_change"],.05186338144743141)
        self.assertEqual(r["s_lat_individual_contrasts_without_mac_range"],[.5057035278775084,1.0])
        a=next(x for x in r["species"] if x["plant"]=="Silene latifolia")
        self.assertEqual((a["n_studies_original"],a["n_studies_strict"]),(5,4))
    def test_dark_V_is_not_same_site_ALAN_measurement(self):
        r=m.compute(self.source,self.original)
        self.assertTrue(r["strict_is_not_independently_measured_local_ALAN_V"])
        self.assertIsNone(r["global_beta_V"])
        self.assertTrue(r["original_species_geographical_model_not_refitted"])
    def test_fails_if_treatment_ratio_altered_without_loss(self):
        rows=[dict(x) for x in self.source]
        one=next(x for x in rows if x["used"]=="True")
        one["open_value"]="0.2"
        with self.assertRaises(AssertionError):m.compute(rows,self.original)
    def test_fails_if_old_study_count_drift(self):
        rows=[dict(x) for x in self.source]
        rows.pop()
        with self.assertRaises(AssertionError):m.compute(rows,self.original)

if __name__=="__main__":
    unittest.main()
