import io
import hashlib
from unittest.mock import patch
import csv
import json
import unittest
import zipfile

from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import (
    parse_zip,PUBLISHED_PAPERS,PUBLISHED_EFFECTS,zenodo_original_author_mirror,ZENODO_EXPECTED_ORIGINAL_MD5
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

    def test_original_zenodo_mirror_validates_real_upstream_checksum(self):
        rawcsv=b"Study,Effect size\nStudyA,0.4\n"
        rawr=b"# genuine-format original source code surrogate for deterministic test\n"
        names={"Sanders_et_al_datafile.csv":rawcsv,"Sanders_meta__rcode.r":rawr}
        md5={name:hashlib.md5(value).hexdigest() for name,value in names.items()}
        record={"metadata":{"related_identifiers":[{"identifier":"10.5061/dryad.wpzgmsbjn"}]},
                "files":[{"key":name,"checksum":"md5:"+md5[name],
                          "links":{"self":"https://zenodo.org/api/records/3998338/files/"+name+"/content"}}
                         for name in names]}
        def fetch(url):
            if url.endswith("/3998338"):
                return json.dumps(record).encode()
            for name,byte in names.items():
                if "/"+name+"/content" in url:
                    return byte
            raise AssertionError(url)
        with patch.dict(ZENODO_EXPECTED_ORIGINAL_MD5,md5,clear=True):
            zip_bytes,files=zenodo_original_author_mirror(fetch)
        self.assertEqual(len(files),2)
        x=parse_zip(zip_bytes,strict=True)
        self.assertEqual(x["verified_original_source_tabular_rows"],1)
        self.assertEqual(x["source_status"],"AUTHENTIC_ORIGINAL_HEDGES_D_CSV_AND_R_CODE_VERIFIED")
        self.assertNotIn("StudyA",str(x))
        record["files"][0]["checksum"]="md5:incorrect"
        with patch.dict(ZENODO_EXPECTED_ORIGINAL_MD5,md5,clear=True):
            with self.assertRaisesRegex(ValueError,"checksum"):
                zenodo_original_author_mirror(fetch)

    def test_original_cr_only_line_endings_are_valid_source_not_fake_csv(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("Sanders_et_al_datafile.csv",
                      b"Study,Effect size,Outcome category\rA,0.35,abundance\rB,-0.12,phenology\r")
            z.writestr("Sanders_meta__rcode.r",b"# test code only\r")
        result=parse_zip(b.getvalue(),strict=True)
        self.assertEqual(result["verified_original_source_tabular_rows"],2)
        self.assertEqual(result["source_original_column_count"],3)
        self.assertNotIn("A,0.35",str(result))

    def test_dont_treat_arbitrary_string_as_original_archival_zip(self):
        with self.assertRaisesRegex(ValueError,"ZIP"):
            parse_zip(b"not-a-real-original-Dryad-file",strict=True)

if __name__=="__main__":
    unittest.main()
