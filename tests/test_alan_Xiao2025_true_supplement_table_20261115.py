import importlib.util, io, zipfile
from pathlib import Path
import unittest
SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_Xiao2025_true_supplement_table_20261115.py"
spec=importlib.util.spec_from_file_location("alan_xiao_table",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class TestXiaoTable(unittest.TestCase):
 def test_original_publisher_hash_pinned(self):
    self.assertTrue(m.URL.endswith("S0006320725005099-mmc1.docx"))
    self.assertEqual(len(m.SHA256),64)
 def test_fake_empty_body_reject(self):
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        z.writestr("word/document.xml",'<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"></w:document>')
    with self.assertRaises(ValueError):m.extract(b.getvalue())
 def test_factor_rows_extracted_without_claim_of_causation(self):
    body='<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:p><w:r><w:t>Table S1. ALAN reproduction</w:t></w:r></w:p><w:tbl><w:tr><w:tc><w:p><w:r><w:t>Factor</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>Effect</w:t></w:r></w:p></w:tc></w:tr><w:tr><w:tc><w:p><w:r><w:t>Artificial light</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>-0.5 p = 0.1</w:t></w:r></w:p></w:tc></w:tr></w:tbl></w:body></w:document>'
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:z.writestr("word/document.xml",body)
    d=m.extract(b.getvalue())
    self.assertEqual(d["document_tables"],1)
    self.assertEqual(d["table_shapes_and_bounded_source_numerics"][0]["row_count"],2)
    self.assertEqual(d["table_shapes_and_bounded_source_numerics"][0]["n_factor_rows_selected"],1)
    self.assertFalse(d["source_published_table_factor_beta_SE_reconstructed"])
if __name__=="__main__":unittest.main()
