"""Minimum source-code dependency audit safeguards."""
import unittest
from scripts.analysis.alan_Macgregor2019_original_script_data_dependency_20261114 import (
    NEEDED, EXPECTED_SCRIPT_GIT_BLOB_SHA1, SCRIPT, AUTHOR
)
class MacgregorEvidence(unittest.TestCase):
    def test_expected_author_source_files_not_imputed(self):
        self.assertEqual(NEEDED,{"SeedSetBinom.csv","SeedSetYN.csv","SeedCount.csv"})
    def test_real_script_sha_and_repo_identity_pinned(self):
        self.assertEqual(EXPECTED_SCRIPT_GIT_BLOB_SHA1,"04136ea82d6c1bbef22f4fd60c41de9907b00550")
        self.assertEqual(AUTHOR,"CallumJMacgregor/BioassayPaper")
        self.assertEqual(SCRIPT,"FinalScripts.R")
if __name__=="__main__":unittest.main()
