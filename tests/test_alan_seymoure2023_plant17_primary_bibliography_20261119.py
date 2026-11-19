import hashlib
import io
import unittest
from unittest.mock import patch
from openpyxl import Workbook
from scripts.acquisition import alan_seymoure2023_plant17_primary_bibliography_20261119 as a

def rec(row,typ,paper,doi):
    return {
        "_original_excel_row":row,
        "interaction_type":typ,"method":"",
        "Authors":paper,"Article Title":"Primary "+paper,
        "Publication Year":"2020","DOI":doi,"DOI Link":"",
        "interaction strength (IS) reduction":"0","IS enhancement":"0",
        "strength sum":"0","temporal mismatch induced":"0",
        "spatial mismatch induced":"0","effect sum":"0",
    }

def test_source():
    rows=[
        rec(2,"pollination","A","10.5555/a"),
        rec(3,"herbivory","A","https://doi.org/10.5555/a"),
        rec(4,"pollination","B","10.5555/b"),
        rec(5,"predator-prey","C","10.5555/c"),
    ]
    footer=rec(6,"","","")
    footer.update({
        "interaction_type":"","Authors":"","Article Title":"",
        "Publication Year":"","DOI":"",
        "interaction strength (IS) reduction":"29",
        "IS enhancement":"28","strength sum":"57",
        "temporal mismatch induced":"7","spatial mismatch induced":"27",
        "effect sum":""
    })
    return rows+[footer]

class Plant17OriginalBibliographyUnitTests(unittest.TestCase):
    def test_shared_study_two_interactions_once(self):
        x=a.source_bibliography(test_source(),strict=False)
        self.assertEqual(len(x),2)
        self.assertEqual(x[0]["doi"],"10.5555/a")
        self.assertEqual(x[0]["original_2023_review_source_interaction_count"],2)
        self.assertEqual(x[0]["original_2023_review_source_types"],["herbivory","pollination"])
        self.assertEqual(x[1]["original_2023_review_source_types"],["pollination"])
        self.assertTrue(all(r["verified_primary_endpoint_effect_size"] is None for r in x))
        self.assertTrue(all(r["verified_primary_endpoint_SE"] is None for r in x))

    def test_original_articles_need_real_doi(self):
        x=test_source()
        x[2]["DOI"]=""
        with self.assertRaisesRegex(ValueError,"lacks DOI"):
            a.source_bibliography(x,strict=False)

    def test_same_study_must_not_have_conflicting_doi(self):
        x=test_source()
        # Same DOI must not silently absorb two conflicting primary titles.
        x[1]["Article Title"]="Different primary paper title"
        with self.assertRaisesRegex(ValueError,"title inconsistency"):
            a.source_bibliography(x,strict=False)

if __name__=="__main__":
    unittest.main()
