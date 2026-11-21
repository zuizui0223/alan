import unittest
from urllib.parse import unquote, urlparse, parse_qs

from scripts.acquisition.alan_seymoure2023_pollination_herbivory_primary_source_access_20261119 import discover


def row(idx,typ,paper,doi,method=""):
    return {
        "_original_excel_row":idx,
        "interaction_type":typ,"method":method,
        "Authors":paper,"Article Title":"Title "+paper,
        "Publication Year":"2021","DOI":doi,"DOI Link":"",
        "interaction strength (IS) reduction":"0",
        "IS enhancement":"0","strength sum":"0",
        "temporal mismatch induced":"0",
        "spatial mismatch induced":"0","effect sum":"0",
    }


def example():
    original=[
        row(2,"pollination","A","10.1234/a","field experiment"),
        row(3,"herbivory","A","10.1234/a"),
        row(4,"herbivory","B",""),
        row(5,"predator-prey","C","10.1234/c"),
    ]
    footer=row(6,"","","")
    footer.update({
        "interaction_type":"",
        "Authors":"","Article Title":"","Publication Year":"",
        "DOI":"",
        "interaction strength (IS) reduction":"29",
        "IS enhancement":"28","strength sum":"57",
        "temporal mismatch induced":"7","spatial mismatch induced":"27",
        "effect sum":"",
    })
    return original+[footer]


def mock_query(url):
    if "/works/" in url:
        doi=unquote(url.split("/works/")[1])
        return {"message":{"DOI":doi,"link":[]}}
    if "europepmc" in url:
        doi=parse_qs(urlparse(url).query)["query"][0][4:]
        return {"resultList":{"result":[
            {"doi":doi,"isOpenAccess":"Y","pmcid":"PMC123"}
        ]}}
    raise AssertionError("Unexpected fabricated metadata URL")


class PlantSourceAcquisitionTests(unittest.TestCase):
    def test_all_source_interaction_types_not_cherry_picked_for_sign(self):
        result=discover(example(),query=mock_query,strict=False)
        self.assertEqual(result["original_priority_interaction_rows"],3)
        self.assertEqual(result["unique_original_publication_keys"],2)
        self.assertEqual(result["original_source_DOI_present_publications"],1)
        self.assertEqual(result["original_method_documented_priority_publications"],1)
        self.assertEqual(result["metadata_access_status"]["europepmc_DOI_exact_OA_reported"],1)
        self.assertFalse("10.1234" in str(result))
        self.assertTrue(result["scientific_stop"]["no_auto_classification_of_randomized_light_experiment"])

    def test_unindexed_doi_is_not_fulltext_unavailable(self):
        def blocked(url):
            if "crossref.org" in url: return {"message":{"DOI":"10.1234/a"}}
            if "europepmc" in url: return {"resultList":{"result":[]}}
            raise AssertionError(url)
        result=discover(example(),query=blocked,strict=False)
        self.assertEqual(result["metadata_access_status"]["europepmc_DOI_NOT_INDEXED"],1)
        self.assertTrue(result["scientific_stop"]["EuropePMC_DOI_not_indexed_is_not_evidence_fulltext_unavailable"])


if __name__=="__main__":
    unittest.main()
