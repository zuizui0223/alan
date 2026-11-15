import importlib.util
import io
import zipfile
import unittest
from pathlib import Path

PATH=Path(__file__).resolve().parents[1] / "scripts/analysis/alan_A_lu2026_source_probe_20261115.py"
spec=importlib.util.spec_from_file_location("source_probe", PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestLuSourceProbe(unittest.TestCase):
    def test_endpoint_is_publisher_or_public_index(self):
        self.assertEqual(len(m.URLS), 8)
        self.assertTrue(all(u.startswith("https://") for u in m.URLS.values()))
        self.assertTrue(all("rtag160" in u for k,u in m.URLS.items() if "metadata" not in k))
    def test_docx_fake_zip_rejected(self):
        self.assertFalse(m.inspect_docx(b"bad")["is_real_docx"])
    def test_actual_docx_structure(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("word/document.xml",
                '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'
                '<w:tbl><w:tr><w:tc><w:p><w:r><w:t>seed set 42 light</w:t></w:r></w:p></w:tc></w:tr></w:tbl>'
                '</w:body></w:document>')
        r=m.inspect_docx(buf.getvalue())
        self.assertTrue(r["is_real_docx"])
        self.assertEqual(r["document_tables"],1)
        self.assertEqual(r["numeric_cells"],1)

if __name__=="__main__":
    unittest.main()
