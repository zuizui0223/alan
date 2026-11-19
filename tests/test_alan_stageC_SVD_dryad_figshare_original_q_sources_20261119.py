import io
import zipfile
import unittest
from scripts.acquisition.alan_stageC_SVD_dryad_figshare_original_q_sources_20261119 import (
  parse_table_header, DRYAD_DOI, FIGSHARE_ID,
)


class OriginalSVDMothFlySourceTests(unittest.TestCase):
    def test_first_original_source_has_bee_fly_not_moth(self):
        self.assertEqual(DRYAD_DOI,"10.5061/dryad.05qfttf1n")
        self.assertEqual(FIGSHARE_ID,26502400)

    def test_original_svd_csv_header_and_no_single_visit_body_records_exported(self):
        text=b"Plant,Guild,Mean SVD,SD,N\nTest flower,Fly,30,4,7\nTest flower,Bee,45,6,8\n"
        d=parse_table_header("Metadata.csv",text)
        self.assertEqual(d["format"],".csv")
        self.assertEqual(d["original_nonempty_data_rows"],2)
        self.assertEqual(d["original_column_names"][0],"Plant")
        self.assertTrue(d["original_biological_values_not_exported"])
        self.assertNotIn("Test flower",str(d))
        self.assertNotIn(",30,",str(d))

    def test_original_zip_source_extracts_no_individual_pollen_grain_rows(self):
        z=io.BytesIO()
        with zipfile.ZipFile(z,"w") as archive:
            archive.writestr("Metadata.csv","Species,Insect,Stigma pollen\nPlantA,Hoverfly,31\n")
        d=parse_table_header("original.zip",z.getvalue())
        self.assertEqual(d["format"],".zip")
        self.assertEqual(len(d["source_members"]),1)
        self.assertEqual(d["source_members"][0]["member_schema"]["original_nonempty_data_rows"],1)
        self.assertNotIn("PlantA",str(d))

    def test_non_source_zip_rejected(self):
        with self.assertRaisesRegex(ValueError,"Not original zip"):
            parse_table_header("original.zip",b"not author actual source")

if __name__=="__main__":
    unittest.main()
