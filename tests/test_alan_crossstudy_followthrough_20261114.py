import math
import unittest
from scripts.analysis.alan_ceiba2022_common_record_associations_20261114 import ranks, spearman, source_ratio_check, association
from scripts.analysis.alan_giavi2020_reported_contrasts_20261114 import number, sections, analyze

class FollowthroughTests(unittest.TestCase):
    def test_tied_ranks(self):
        self.assertEqual(ranks([3,1,1,5]),[3,1.5,1.5,4])
    def test_signed_rank_relationship(self):
        self.assertAlmostEqual(spearman([1,1,3,4],[2,2,9,12]),1.0)
        self.assertAlmostEqual(spearman([1,1,3,4],[-2,-2,-9,-12]),-1.0)
    def test_no_unequal_samples_or_constant_outcomes(self):
        with self.assertRaises(ValueError):spearman([1,1,1],[1,2,3])
        with self.assertRaises(ValueError):spearman([1,2],[1,2,3])
    def test_ratio_check_not_total_seed_invention(self):
        rows=[{'SeedsN':'10','Ovules':'20','SeedSet':'.5'},
              {'SeedsN':'3','Ovules':'9','SeedSet':'.333'}]
        r=source_ratio_check(rows,'SeedsN','Ovules','SeedSet')
        self.assertEqual(r['within_0p001'],2)
        self.assertTrue(r['identity_does_not_prove_shared_flower_sampling'])
    def test_leaveout_not_independent_replication(self):
        r=association([1,2,3,4,5],[2,4,6,8,10])
        self.assertEqual(r['n'],5)
        self.assertTrue(r['p_and_cluster_CI_not_estimated'])
        self.assertTrue(r['leave_one_out_not_independent_replications'])
    def test_source_unicode_signed_numbers(self):
        self.assertAlmostEqual(number('− 1.90'),-1.9)
    def test_original_reference_orientation_and_denominator_caution(self):
        rows=[['(a)'],['Dark (intercept)','0.43','.656','.7','.512'],
              ['(b)'],['Adjacent (intercept)','.24','.536','.4','.659'],
              ['Illuminated','− 1.14','.363','− 3.1','.002'],
              ['Dark','− 1.90','.576','− 3.3','< 0.001'],
              ['(c)'],['Adjacent (intercept)','− .74','.464','− 1.6','.109'],
              ['Dark','1.04','.499','2.1','.037']]
        source={'original_published_model_tables':[{'label':'Table 1','rows':rows}],
                'article_doi':'10.1038/s41598-020-68667-y','article_xml_sha256':'test'}
        r=analyze(source)
        self.assertAlmostEqual(r['published_predation_contrasts'][0]['odds_ratio'],math.exp(1.9))
        self.assertEqual(r['published_predation_contrasts'][0]['contrast'],'Adjacent / Dark')
        self.assertFalse(r['published_predation_contrasts'][0]['new_raw_data_model'])
        self.assertTrue(r['relative_reproduction_published_coefficient']['denominator_semantics_need_author_code'])
    def test_incomplete_source_table_rejected(self):
        with self.assertRaises(ValueError):sections([['(a)']])

if __name__=='__main__':unittest.main()
