import importlib.util
import unittest
from pathlib import Path

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_P1_original_smallfruit_diel_dependence_20261115.py"
s=importlib.util.spec_from_file_location("p1crops",SRC)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

class TestP1SmallFruit(unittest.TestCase):
    def test_original_hashes_pinned(self):
        self.assertEqual(len(m.FILES),2)
        self.assertEqual(m.FILES["fruit_set"]["expected_rows"],80)
        self.assertEqual(m.FILES["fruit_quality"]["expected_rows"],223)
        self.assertEqual(len(m.FILES["fruit_set"]["sha256"]),64)
    def test_arms_mapped_without_assumption_of_guild(self):
        self.assertEqual(m.normalize_treatment("Day + Night"),"day_night_open")
        self.assertEqual(m.normalize_treatment("Night"),"night_only")
        self.assertEqual(m.normalize_treatment("Closed"),"closed")
        self.assertEqual(m.normalize_treatment("some_unknown"),"UNRESOLVED(some_unknown)")
    def test_nocturnal_marginal_positive_and_negative(self):
        # One row per matched treatment with 10 initial flowers each.
        rows=[
          {"Plant_sp":"Red","Year":"2020","Treatment":"Day","Fruit_set_success":"4","Fruit_set_failure":"6"},
          {"Plant_sp":"Red","Year":"2020","Treatment":"Day + Night","Fruit_set_success":"8","Fruit_set_failure":"2"},
          {"Plant_sp":"Black","Year":"2024","Treatment":"Day","Fruit_set_success":"6","Fruit_set_failure":"4"},
          {"Plant_sp":"Black","Year":"2024","Treatment":"Day + Night","Fruit_set_success":"3","Fruit_set_failure":"7"}
        ]
        d=m.reduce_source(rows,[])
        index={(x["taxon_original"],x["year"]):x["V_fruit_night_access"] for x in d["new_night_access_V_by_taxon_year"]}
        self.assertAlmostEqual(index[("Red","2020")],.5)
        self.assertAlmostEqual(index[("Black","2024")],-1)
        self.assertFalse(d["ALAN_risk_ratio_estimated"])
    def test_real_zeros_are_allowed_missing_not_imputed(self):
        rows=[
          {"Plant_sp":"R","Year":"2020","Treatment":"Day","Fruit_set_success":"0","Fruit_set_failure":"2"},
          {"Plant_sp":"R","Year":"2020","Treatment":"Day + Night","Fruit_set_success":"1","Fruit_set_failure":"1"}
        ]
        out=m.reduce_source(rows,[])
        self.assertEqual(out["new_night_access_V_by_taxon_year"][0]["V_fruit_night_access"],1.0)
        self.assertIsNone(out["new_night_access_V_by_taxon_year"][0]["SE_or_CI"])
    def test_reject_noninteger_source_count(self):
        with self.assertRaises(ValueError):
            m.int_count("NaN")
        with self.assertRaises(ValueError):
            m.int_count("2.4")

if __name__=="__main__":
    unittest.main()
