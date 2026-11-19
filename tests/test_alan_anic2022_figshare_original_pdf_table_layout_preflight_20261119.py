import io
import unittest
from unittest.mock import patch

from scripts.acquisition.alan_anic2022_figshare_original_pdf_table_layout_preflight_20261119 import (
    masking, study_pdf_preflight,
)


class Anic2022OriginalNumericPdfSourceTests(unittest.TestCase):
    def test_numeric_data_and_plot_code_are_redacted_from_headers(self):
        raw="Trough P5 White LED 2016 Lotus pedunculatus -13.42 0.362"
        txt=masking(raw)
        self.assertIn("White LED",txt)
        self.assertIn("Lotus pedunculatus",txt)
        self.assertNotIn("P5",txt)
        self.assertNotIn("2016",txt)
        self.assertNotIn("-13.42",txt)
        self.assertNotIn("0.362",txt)

    def test_dryad_original_doi_not_confused_with_figshare_plant_pdf(self):
        from scripts.acquisition.alan_anic2022_figshare_original_pdf_table_layout_preflight_20261119 import DOI, FILE_ID, FILE_SHA256
        self.assertEqual(DOI,"10.1002/ece3.9157")
        self.assertEqual(FILE_ID,36265164)
        self.assertEqual(FILE_SHA256,"baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0")

    def test_reject_fake_PDF_and_wrong_author_original_SHA(self):
        with self.assertRaisesRegex(ValueError,"SHA mismatch"):
            study_pdf_preflight(b"%PDF- fake",strict=True)

if __name__=="__main__":
    unittest.main()
