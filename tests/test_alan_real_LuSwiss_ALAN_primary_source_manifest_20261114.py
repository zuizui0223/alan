#!/usr/bin/env python3
"""ALAN source-manifest reader guards, never dummy field experiment effects."""
import sys,unittest,io,zipfile,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_real_LuSwiss_ALAN_primary_source_manifest_20261114 import (
 manifest_zip_bytes,PUBLISHED_SUPP,SWISS_ENDPOINTS
)

class OriginalPrimarySourceManifest(unittest.TestCase):
    def test_only_valid_original_uri_not_invented(self):
        self.assertIn("PMC13645961",PUBLISHED_SUPP[0])
        self.assertIn("xjmn5jhxyn",SWISS_ENDPOINTS[0])
        self.assertTrue(all(x.startswith("https://") for x in PUBLISHED_SUPP+SWISS_ENDPOINTS))

    def test_real_zip_header_schema_with_no_rows(self):
        f=io.BytesIO()
        with zipfile.ZipFile(f,"w") as z:
            z.writestr("Supplementary_Readme.txt","title and methods")
        res=manifest_zip_bytes(f.getvalue(),"synthetic")
        self.assertEqual(res["status"],"valid_zip")
        self.assertEqual(len(res["original_entries"]),1)
        self.assertNotIn("rows",res)

    def test_nonzip_does_not_claim_success(self):
        res=manifest_zip_bytes(b'no empirical field data',"test")
        self.assertEqual(res["status"],"not_zip")

if __name__=="__main__":
    unittest.main()
