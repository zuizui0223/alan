"""Carannante four spatial-pair source discovery must not claim sampled data."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"scripts/acquisition/alan_carannante2021_boris_source_gate_20261124.py"
spec=importlib.util.spec_from_file_location("boris_source",FILE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestBernSource(unittest.TestCase):
    def test_access_probes_allowlisted_and_no_data_proven_without_bytes(self):
        self.assertEqual(len(m.SOURCES),4)
        self.assertTrue(all(x.startswith("https://boris") for x in m.SOURCES))
        self.assertEqual(m.DOI,"10.1002/2688-8319.12053")
        self.assertTrue(all(x.startswith("https://") for x in m.SOURCES))
        with self.assertRaises(ValueError):
            m.probe("http://evil.test/carannante")
        with self.assertRaises(ValueError):
            m.probe("https://github.com/zuizui0223/alan")
    def test_public_dataset_name_hints_are_not_original_dataset(self):
        obj={"_embedded":{"searchResult":{"_embedded":{"objects":[
            {"_embedded":{"indexableObject":{"name":"Carannante et al adult aquatic insects 2021"}}},
            {"name":"A different document"},
        ]}}}}
        self.assertEqual(m.extract_name_hints(obj),
                         ["Carannante et al adult aquatic insects 2021"])
    def test_html_extraction_is_only_filename_suffix_gate(self):
        parser=m.CandidateRefs()
        parser.feed('<a href="/download/original.xlsx">Excel</a><a href="/login">Login</a>')
        self.assertEqual(parser.urls,["/download/original.xlsx"])

if __name__=="__main__":
    unittest.main()
