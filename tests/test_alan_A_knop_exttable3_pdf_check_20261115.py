import importlib.util
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_A_knop_exttable3_pdf_check_20261115.py"
spec=importlib.util.spec_from_file_location("knop_pdf",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class KnopTableTest(unittest.TestCase):
    def test_source_bound(self):
        self.assertEqual(m.LIMIT,125_000_000)
        self.assertIn("planesynormas.mma.gob.cl",m.URL)
    def test_public_source_only(self):
        self.assertTrue(m.URL.endswith(".pdf"))
    def test_extract_empty_pdf_reject(self):
        with self.assertRaises(Exception):
            m.extract_table_text(b"abc")
if __name__=="__main__":
    unittest.main()
