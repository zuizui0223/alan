#!/usr/bin/env python3
"""Dryad file source schema and DOI pin only; original API data separate."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_moth23_Dryad_original_source_manifest_20261114 import (
    DOI,ENC,BASE,hrefs,objects,safe_file
)
class DryadManifest(unittest.TestCase):
 def test_versioned_doi_and_encoded_slash(self):
  self.assertEqual(DOI,"doi:10.5061/dryad.jq2bvq8ps")
  self.assertIn("%2F",ENC)
  self.assertEqual(BASE,"https://datadryad.org/api/v2")
 def test_original_file_metadata_no_video_rows(self):
  x={"id":7,"path":"Activity_MainData.csv","size":77,"digest":"a","digestType":"sha256",
     "Species":"sensitive","Moth_ID":"private"}
  y=safe_file(x)
  self.assertEqual(y["path"],"Activity_MainData.csv")
  self.assertNotIn("Species",y)
  self.assertNotIn("Moth_ID",y)
 def test_links(self):
  z={"_links":{"file":{"href":"https://datadryad.org/api/v2/versions/7/files"}}}
  self.assertEqual(list(hrefs(z)),["https://datadryad.org/api/v2/versions/7/files"])
  self.assertTrue(any(v.get("href") for v in objects(z)))
if __name__=="__main__":unittest.main()
