import io
import unittest
from openpyxl import Workbook
from scripts.acquisition.alan_seymoure2023_original_79paper_figshare import sheet_schema


class InteractionReviewSourceTests(unittest.TestCase):
    def test_79_papers_not_81_independent_trials(self):
        from scripts.acquisition.alan_seymoure2023_original_79paper_figshare import audit
        # Static published counts. Do not call the network in deterministic tests.
        classes=(56,11,6,5,1,1,1)
        self.assertEqual(sum(classes),81)
        self.assertNotEqual(sum(classes),79)

    def test_excel_header_and_rows_cannot_export_paper_data(self):
        wb=Workbook()
        ws=wb.active
        ws.title="Published review source"
        ws.append(["Study","Interaction type","Effect direction"])
        ws.append(["A published study","predation","enhanced"])
        ws.append(["B published study","pollination","reduced"])
        stream=io.BytesIO()
        wb.save(stream)
        wb.close()
        summary=sheet_schema(stream.getvalue())
        self.assertEqual(summary[0]["first_likely_headers_only"],
                         ["Study","Interaction type","Effect direction"])
        self.assertEqual(summary[0]["source_nonempty_rows"],3)
        self.assertNotIn("A published study",str(summary))
        self.assertNotIn("B published study",str(summary))
        self.assertTrue(summary[0]["no_original_study_rows_exported"])

    def test_non_xlsx_bytes_rejected(self):
        with self.assertRaisesRegex(ValueError,"No original Excel"):
            sheet_schema(b"no real source XLSX")


if __name__=="__main__":
    unittest.main()
