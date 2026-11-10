"""Only synthetic calendar and parsing guards for the true GloPL × historical PCNL join."""
import unittest
from scripts.analysis import alan_glopl_prior_PCNL_6year_shard_20261110 as a


class GloPLPcnlGuards(unittest.TestCase):
    def test_1992_to_2014_are_23_disjoint_original_annual_pcNL_rasters(self):
        got=sum([a.years(s) for s in range(6)],[])
        self.assertEqual(sorted(got),list(range(1992,2015)))
        self.assertEqual(len(set(got)),23)
    def test_wrong_shard_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"six frozen"):
            a.years(7)
    def test_preobservational_source_checksum_is_complete(self):
        self.assertEqual(len(a.SOURCE_SHA),64)
        self.assertEqual(a.SOURCE_ROWS,1419)
        self.assertIn("experiment",a.SCHEMA)
if __name__=="__main__": unittest.main()
