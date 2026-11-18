import io
import unittest
from openpyxl import Workbook
from scripts.acquisition.alan_eckhartt2022_original_dryad_7lamps_20261118 import workbook_schema


class EckharttOriginalDryadSourceGateTests(unittest.TestCase):
    def workbook(self):
        w=Workbook()
        s=w.active
        s.title="Original Data"
        s.append(["Lamppost","TimeCondition","Distance_m","Predated"])
        s.append(["lamp_A","before",1,1])
        s.append(["lamp_A","after",1,0])
        s.append(["lamp_B","after",2,1])
        b=io.BytesIO()
        w.save(b)
        w.close()
        return b.getvalue()

    def test_detect_only_original_id_counts_not_lamp_ids(self):
        t=workbook_schema(self.workbook())
        self.assertEqual(len(t),1)
        self.assertEqual(t[0]["source_rows_after_header"],3)
        self.assertEqual(t[0]["header"][0],"Lamppost")
        self.assertEqual(t[0]["possible_lamppost_ID_unique_counts_only"]["0"],2)
        self.assertNotIn("lamp_A",str(t))
        self.assertNotIn("lamp_B",str(t))

    def test_no_fabricated_non_xlsx_data(self):
        with self.assertRaisesRegex(ValueError,"not an XLSX"):
            workbook_schema(b"this is a placeholder, not original data")


if __name__=="__main__":
    unittest.main()
