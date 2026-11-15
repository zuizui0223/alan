import unittest
from scripts.analysis.alan_Grognuz2026_slug_2x2_public_raw_probe_20261115 import DATASET_ID,DOI,ENDPOINTS,request
class PublicSourceProbeContracts(unittest.TestCase):
    def test_published_dataset_exact_id_and_version(self):
        self.assertEqual(DATASET_ID,"drhzznh4ps")
        self.assertEqual(DOI,"10.17632/drhzznh4ps.1")
        self.assertTrue(any("version=1" in v for v in ENDPOINTS.values()))
    def test_official_api_and_public_files_not_private_auth(self):
        self.assertTrue(ENDPOINTS["public_dataset_v1"].startswith("https://api.data.mendeley.com/datasets/"))
        self.assertTrue(ENDPOINTS["public_dataset_files"].startswith("https://api.data.mendeley.com/datasets/publics/"))
    def test_bounded_receipt_on_bad_url(self):
        x,meta=request("invalid-url",limit=100)
        self.assertIsNone(x)
        self.assertEqual(meta["status"],"URL_ERROR")
        self.assertEqual(meta["error_kind"],"ValueError")
if __name__=="__main__":unittest.main()
