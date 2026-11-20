"""Unit and scientific accuracy tests of coordinate-source coverage."""
import importlib.util
import pathlib
import sys
import unittest
import pandas as pd
S=pathlib.Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(S))
spec=importlib.util.spec_from_file_location("recover",S/"alan_recover_frozen187_pair_coordinates_20261120.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class TestRecover(unittest.TestCase):
    def test_no_conflicting_site_ids(self):
        z=pd.DataFrame([{"unit_id":1,"longitude":-1.,"latitude":50.},
                        {"unit_id":1,"longitude":-2.,"latitude":50.}])
        with self.assertRaisesRegex(ValueError,"CONFLICT"):
            m.validate_coordinates(z)
    def test_partial_coverage_never_counted_as_complete(self):
        p=pd.DataFrame([{"bright_site_id":1,"comparison_site_id":2,"block":1},
                        {"bright_site_id":3,"comparison_site_id":4,"block":2}])
        coords=pd.DataFrame([{"unit_id":1},{"unit_id":2},{"unit_id":3}])
        a=m.pair_coverage(p,coords)
        self.assertEqual(a["matched_pairs_both_coordinates"],1)
        self.assertEqual(a["unique_original_matched_sites"],4)
    def test_real_public_site_sources_sha_and_coordinate_consistency(self):
        c=m.known_repository_coordinates()
        self.assertGreaterEqual(len(c),18)
        self.assertFalse(c.unit_id.duplicated().any())
        self.assertTrue(c.source_count.ge(1).all())
    def test_original_matched_pair_count(self):
        p,_=m.matched_ids()
        self.assertEqual(len(p),187)
        self.assertEqual(p.block.nunique(),30)
if __name__=="__main__":
    unittest.main()
