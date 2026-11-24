"""Source-authentic, independent source row and design-fidelity tests."""
import math, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/analysis'))
from alan_sevenponds_family_sensitivity_20261123 import run,source,interaction, violet_moderator

class SevenPondsOriginalSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=run()
    def test_exact_30_lamp_10_cells(self):
        self.assertEqual(self.data['true_independent_lighting_units'],30)
        self.assertEqual(self.data['original_taxa_rows'],193)
        self.assertEqual(self.data['families'],35)
        self.assertEqual(self.data['original_invertebrates'],4699)
    def test_published_original_anova_reproduced(self):
        self.assertAlmostEqual(self.data['F4_20_profiles']['full_raw']['F'],9.680484236634571,places=8)
        self.assertAlmostEqual(self.data['F4_20_profiles']['full_raw']['p_nominal'],0.00015950472855757,places=9)
    def test_taxon_identity_and_abundance(self):
        top=self.data['top_two_families']
        self.assertEqual([v['name'] for v in top],['caenidae','leptoceridae'])
        self.assertEqual([v['count'] for v in top],[2860.,686.])
        self.assertAlmostEqual(self.data['top_two_fraction'],3546/4699,places=10)
    def test_transform_and_taxon_removal_change_interaction(self):
        s=self.data['F4_20_profiles']
        self.assertAlmostEqual(s['full_log1p']['F'],6.757759234236177,places=9)
        self.assertAlmostEqual(s['top2_excluded_raw']['F'],4.2930071984721625,places=9)
        self.assertAlmostEqual(s['top2_excluded_log1p']['F'],2.760511979724383,places=9)
        self.assertGreater(s['top2_excluded_log1p']['p_nominal'],0.05)
    def test_exhaustive_two_family_scan_no_filtered_suppression(self):
        self.assertEqual(self.data['pairs_total'],595)
        self.assertEqual(self.data['pairs_p_at_least_0_05'],2)
        self.assertEqual(self.data['pairs_nonsignificant'],[['cicadellidae','leptoceridae'],['caenidae','leptoceridae']])
    def test_violet_littoral_riparian_ratio(self):
        self.assertAlmostEqual(self.data['violet_full']['ratio_of_RRs'],10.748905565978736,places=10)
        self.assertAlmostEqual(self.data['violet_excluding_top_two']['ratio_of_RRs'],5.2435064935064934,places=10)
    def test_carannante_exact_support_only_idealised(self):
        z=self.data['carannante_design_idealised_small_n']
        self.assertEqual(z['spatial_blocks'],2)
        self.assertEqual(z['two_sided_min_p_if_two_blocks_exchangeable'],0.5)
        self.assertAlmostEqual(z['two_sided_min_p_if_four_pairs_fully_exchangeable'],1/3)
        self.assertIn('Hypothetical',z['note'])

if __name__=='__main__':unittest.main()
