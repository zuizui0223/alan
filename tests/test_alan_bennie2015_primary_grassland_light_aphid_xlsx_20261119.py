import io
import unittest
from openpyxl import Workbook
from scripts.acquisition.alan_bennie2015_primary_grassland_light_aphid_xlsx_20261119 import workbook_schema

class SourcePrimaryBennie2015Tests(unittest.TestCase):
    def sample(self):
        wb=Workbook()
        ws=wb.active
        ws.title="original tab"
        ws.append(["plot ID","lighting","aphids","flower heads"])
        ws.append(["PLOT_SENSITIVE_ABC","LED",40,3])
        ws.append(["PLOT_SENSITIVE_DEF","DARK",20,6])
        b=io.BytesIO()
        wb.save(b)
        wb.close()
        return b.getvalue()

    def test_no_original_plot_data_export_and_real_shape(self):
        x=workbook_schema(self.sample())
        self.assertEqual(len(x),1)
        self.assertEqual(x[0]["data_nonempty_rows"],3)
        self.assertEqual(x[0]["sheet"],"original tab")
        self.assertTrue(x[0]["original_source_plot_IDs_not_exported"])
        self.assertNotIn("PLOT_SENSITIVE_ABC",str(x))
        self.assertNotIn("PLOT_SENSITIVE_DEF",str(x))
        self.assertNotIn("aphids",str(x))

    def test_no_fake_spreadsheet(self):
        with self.assertRaises(ValueError):
            workbook_schema(b"not genuine source")

if __name__=="__main__":
    unittest.main()
