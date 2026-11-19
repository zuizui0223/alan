import io
import csv
import json
import unittest
import zipfile

from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import (
    parse_zip,PUBLISHED_PAPERS,PUBLISHED_EFFECTS
)

class Sanders2020OriginalALANReviewSourceTests(unittest.TestCase):
    def make_source(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("Sanders_et_al_datafile.csv",
                "Study,Effect size,Outcome category,Light type,Trophic guild\n"
                "FictitiousStudyA,0.34,abundance,LED,spider\n"
                "FictitiousStudyA,-0.42,phenology,LED,plant\n"
                "FictitiousStudyB,0.10,reproduction,amber,plant\n")
            z.writestr("Sanders_meta__rcode.r",
                "# synthetic test-only R code, NOT author original")
        return b.getvalue()

    def test_source_csv_and_code_sha_and_no_original_study_rows_export(self):
        x=parse_zip(self.make_source(),strict=True)
        self.assertEqual(x["published_original_papers"],126)
        self.assertEqual(x["published_original_effects"],1304)
        self.assertEqual(x["source_status"],"AUTHENTIC_ORIGINAL_HEDGES_D_CSV_AND_R_CODE_VERIFIED")
        self.assertEqual(x["verified_original_source_tabular_rows"],3)
        self.assertFalse(x["original_CSV_exact_1304_effect_rows_matches_published"])
        self.assertIn("Outcome category",x["source_original_column_names"])
        self.assertEqual(x["author_effect_metric"],"Hedges d, NOT lnRR")
        self.assertNotIn("FictitiousStudyA",str(x))
        self.assertNotIn("FictitiousStudyB",str(x))
        self.assertFalse(x["real_water_wetland_moderator_source_validated"])

    def test_original_data_source_provenance_not_a_meta_ready_lnRR(self):
        self.assertEqual(PUBLISHED_PAPERS,126)
        self.assertEqual(PUBLISHED_EFFECTS,1304)
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("Sanders_et_al_datafile.csv","Study,d\nA,0.5\n")
        with self.assertRaisesRegex(ValueError,"author CSV and code"):
            parse_zip(b.getvalue(),strict=True)

    def test_dont_treat_arbitrary_string_as_original_archival_zip(self):
        with self.assertRaisesRegex(ValueError,"ZIP"):
            parse_zip(b"not-a-real-original-Dryad-file",strict=True)

if __name__=="__main__":
    unittest.main()
