"""Source-language and checksum gates for independent day/night exclusion studies."""
from __future__ import annotations
import hashlib
import tempfile
import unittest
from pathlib import Path
import pandas as pd

from scripts.analysis import alan_diel_exclusion_experimental_source_overlap_20261110 as audit


class TestDielSource(unittest.TestCase):
    def test_original_git_object_hash_form(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"blob.csv"
            raw=b"plant_species,life_form\nSalix caprea,tree\n"
            p.write_bytes(raw)
            expected=hashlib.sha1((f"blob {len(raw)}\0").encode()+raw).hexdigest()
            self.assertEqual(audit.git_blob_sha(p),expected)

    def test_normalized_names_do_not_imply_synonyms_or_pollinator_identity(self):
        x=pd.Series([" Salix   caprea ","Adenophora\tcapillaris","Unknown plant"])
        self.assertEqual(audit.normalize_name(x).tolist(),
                         ["Salix caprea","Adenophora capillaris","Unknown plant"])
        self.assertFalse(hasattr(audit,"is_obligate_moth_pollinated"))

    def test_original_source_is_limited_and_not_a_plant_absence(self):
        self.assertEqual(audit.TRAIT_BLOB,"18381e550a7e405fe4f198909f0bb381bbecf1e0")
        self.assertEqual(audit.STUDY_BLOB,"39a82d966afcdde670603add1c881ace8791ffe9")
        self.assertEqual(audit.KEY,["island","latbin","lonbin","year"])

    def test_missing_file_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/"fake.csv"
            a.write_text("plant_species\nSalix caprea\n")
            with self.assertRaisesRegex(ValueError,"Git blob changed"):
                audit.experimental_records(a,a)


if __name__=="__main__":
    unittest.main()
