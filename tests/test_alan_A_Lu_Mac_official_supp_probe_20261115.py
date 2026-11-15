import io, zipfile, importlib.util
from pathlib import Path
import unittest
path=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_A_Lu_Mac_official_supp_probe_20261115.py"
spec=importlib.util.spec_from_file_location("suppprobe",path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestOfficialSourceProbe(unittest.TestCase):
    def test_source_hosts_official(self):
        self.assertTrue(any("academic.oup.com" in x for x in m.SOURCES.values()))
        self.assertTrue(any("jpe.ac.cn" in x for x in m.SOURCES.values()))
        self.assertTrue(any("esajournals.onlinelibrary.wiley.com" in x for x in m.SOURCES.values()))
    def test_invalid_docx_not_assumed_data(self):
        self.assertFalse(m.docx_summary(b"BAD")["valid_word_docx"])
    def test_valid_docx_with_S4b_detected_without_copying_table(self):
        b=io.BytesIO()
        xml='<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:tbl><w:tr><w:tc><w:p><w:r><w:t>Table S4b seed 12</w:t></w:r></w:p></w:tc></w:tr></w:tbl></w:body></w:document>'
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("word/document.xml",xml)
        d=m.docx_summary(b.getvalue())
        self.assertEqual(d["table_count"],1)
        self.assertTrue(d["has_table_S4b_literal"])
        self.assertTrue(d["has_reproductive_seed_text"])
        self.assertEqual(d["numeric_table_cell_count"],1)
        self.assertFalse(d["raw_plant_rows_recovered"])
if __name__=="__main__":
    unittest.main()
