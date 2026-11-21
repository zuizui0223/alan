"""Source access must not invent mesocosm exposure or numeric biomass."""
import unittest
from scripts.acquisition.alan_xue2026_original_sciencedb_source_gate import (
    scan_meta, source_schema, execute
)

class XueSourceGate(unittest.TestCase):
    def test_only_explicit_original_files(self):
        records={
            "doi":"10.57760/sciencedb.35715",
            "metadata":{"name":"experiment3.csv","fileUrl":"https://www.scidb.cn/real/file.csv","size":1200},
            "abstract":"herbivory native plants changed"
        }
        x=scan_meta(records)
        self.assertEqual(len(x),1)
        self.assertEqual(x[0]["reported_name"],"experiment3.csv")
        self.assertEqual(x[0]["metadata_size"],1200)

    def test_csv_headers_and_unit_are_source_not_guessed(self):
        x=source_schema(b"species,zone,light,biomass\nA,z1,lit,0.5\nB,z2,dark,0.1\n","source.csv")
        self.assertEqual(x["rows"],2)
        self.assertIn("zone",x["columns"])
        self.assertEqual(x["source_factor_levels"]["light"],["dark","lit"])

    def test_html_is_not_parsed_as_original_csv(self):
        z=source_schema(b"<html><body>title only</body></html>","source.csv")
        self.assertEqual(z["type"],"unrecognized_original_CSV")

if __name__=="__main__":
    unittest.main()
