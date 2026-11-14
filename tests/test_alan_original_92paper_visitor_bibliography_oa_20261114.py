#!/usr/bin/env python3
"""Pure safety tests for DOI source normalization and OA metadata-only gate."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_92paper_visitor_bibliography_oa_20261114 import (
    normalize_identifier, metadata_for, year_of_crossref
)


class OriginalPaperOAGate(unittest.TestCase):
    def test_real_original_doi_url_and_literal_are_same_identity(self):
        for s in (
            "https://doi.org/10.1111/boj.12183",
            "10.1111/boj.12183",
            " DOI: 10.1111/boj.12183  ",
            "https://dx.doi.org/10.1111/boj.12183",
        ):
            self.assertEqual(normalize_identifier(s), ("doi_syntax", "10.1111/boj.12183"))

    def test_source_non_doi_is_never_guessed(self):
        self.assertEqual(normalize_identifier(""), ("missing", None))
        self.assertEqual(normalize_identifier("NA"), ("missing", None))
        self.assertEqual(normalize_identifier("NA - have PDF")[0], "invalid_literal")
        self.assertEqual(normalize_identifier("10.590/0001-3765201520140340")[0],
                         "invalid_literal")
        self.assertEqual(normalize_identifier(
            "https://journals.flvc.org/selbyana/article/view/121717")[0],
            "non_doi_url")

    def test_metadata_only_does_not_identify_animal_visitor(self):
        def good(url):
            if "crossref" in url:
                return {"message": {"DOI": "10.1111/boj.12183",
                                    "title": ["Pollination test"],
                                    "issued": {"date-parts": [[2014]]}}}, "ok"
            return {"doi": "https://doi.org/10.1111/boj.12183",
                    "open_access": {"is_oa": True, "oa_status": "green"},
                    "best_oa_location": {"landing_page_url":
                                         "https://example.edu/pub/1"},
                    "has_fulltext": True}, "ok"
        x = metadata_for("10.1111/boj.12183", fetch=good)
        self.assertTrue(x["openalex_oa_record_claim"])
        self.assertTrue(x["candidate_oa_link_NOT_tested_for_article_text"])
        self.assertTrue(
            x["actual_visitor_taxa_or_effective_pollination_NOT_inferred_from_abstract"])
        self.assertEqual(x["crossref_publication_year"], 2014)

    def test_wrong_doi_fails_closed(self):
        def wrong(url):
            if "crossref" in url:
                return {"message": {"DOI": "10.4321/mismatched",
                                    "title": ["Misidentified article"]}}, "ok"
            return {"doi": "https://doi.org/10.4321/mismatched",
                    "open_access": {"is_oa": True},
                    "best_oa_location": {}}, "ok"
        x = metadata_for("10.1111/boj.12183", fetch=wrong)
        self.assertEqual(x["crossref_lookup"], "wrong_doi")
        self.assertIsNone(x["crossref_title"])
        self.assertEqual(x["openalex_lookup"], "wrong_doi")

    def test_year_absence_is_unknown(self):
        self.assertIsNone(year_of_crossref({}))
        self.assertEqual(year_of_crossref({"published": {"date-parts": [[2017, 5]]}}),
                         2017)


if __name__ == "__main__":
    unittest.main()
