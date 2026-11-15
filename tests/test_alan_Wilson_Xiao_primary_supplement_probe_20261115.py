import importlib.util
import unittest
import io
import zipfile
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_Wilson_Xiao_primary_supplement_probe_20261115.py"
spec=importlib.util.spec_from_file_location("alan_xiao_wilson",SOURCE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestProbe(unittest.TestCase):
    def test_known_publications_and_scale(self):
        self.assertEqual(len(m.PAPERS),2)
        self.assertEqual(m.PAPERS["wilson2021"]["doi"],"10.1093/icb/icab010")
        self.assertEqual(m.PAPERS["xiao2025"]["doi"],"10.1016/j.biocon.2025.111472")
        self.assertTrue("40" in m.PAPERS["xiao2025"]["correct_unit"])
    def test_original_repro_scale_not_attributed_to_night_access(self):
        self.assertIn("inflorescence",m.PAPERS["xiao2025"]["plant_outcome"])
        self.assertIn("plant",m.PAPERS["wilson2021"]["plant_outcome"])
    def test_zip_signatures_are_verified(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("table_s1.csv","testvalue")
        entry={"content_type":"application/zip"}
        m.inspect(b.getvalue(),entry)
        self.assertTrue(entry["zip_valid"])
        self.assertTrue(entry["plausible_experimental_csv"])
        self.assertFalse(entry["genuine_reproduction_cell_rows_extracted"])
    def test_html_extracts_only_relevant_urls(self):
        entry={"content_type":"text/html"}
        m.inspect(b'<html><a href="icab010_supplementary_data.zip">x</a><a href="/privacy">privacy</a></html>',entry)
        self.assertIn("icab010_supplementary_data.zip",entry["relevant_public_html_links"])
        self.assertNotIn("/privacy",entry["relevant_public_html_links"])

if __name__=="__main__":
    unittest.main()
