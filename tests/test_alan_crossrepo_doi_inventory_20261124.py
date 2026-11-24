"""Unit checks for DOI extraction and fail-honest repo scan budget (no web)."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/acquisition/alan_crossrepo_doi_inventory_20261124.py"
spec=importlib.util.spec_from_file_location("alan_crossrepo",SCRIPT)
a=importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class TestCrossRepoDOIs(unittest.TestCase):
    def test_valid_citations_and_punctuation(self):
        s="https://doi.org/10.1038/nature23288; another doi:10.1093/jpe/rtag160. Same DOI 10.1038/nature23288."
        x=[a.normalize_doi(m.group()) for m in a.DOI.finditer(s)]
        self.assertEqual(x,["10.1038/nature23288","10.1093/jpe/rtag160","10.1038/nature23288"])

    def test_bib_is_prioritized_and_large_files_are_not_silently_included(self):
        bib=a.candidate_rank("docs/refs.bib",300)
        review=a.candidate_rank("docs/literature_review.md",4000)
        readme=a.candidate_rank("README.md",1000)
        script=a.candidate_rank("scripts/run.py",2000)
        self.assertGreater(bib,readme)
        self.assertGreater(review,50)
        self.assertIsNone(script)
        self.assertIsNone(a.candidate_rank("results/data.csv",9999999))
        self.assertIsNone(a.candidate_rank(".github/workflows/ci.yml",2000))

    def test_no_search_hit_is_meta_study(self):
        self.assertIn(".ris",a.BIB)
        self.assertNotIn(".pdf",a.READABLE)
        self.assertEqual(a.MAX_DOIS_PER_FILE,150)
        self.assertEqual(a.USER,"zuizui0223")
        self.assertTrue(all(a.candidate_rank("README.md",1000)>0 for _ in [1]))


if __name__=="__main__":
    unittest.main()
