import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts' / 'analysis'))
import alan_moth_source_bridge as mod


class RealBridgeLogicTest(unittest.TestCase):
    def test_record_only_is_not_absence(self):
        g = pd.DataFrame({
            'plant_name': ['Flora alba', 'Flora alba', 'Flora rubra'],
            'guild': ['moth_family_visit', 'butterfly_family_visit', 'butterfly_family_visit'],
            'relation': ['visitsFlowersOf', 'visitsFlowersOf', 'pollinates']})
        a = mod.positive_guilds(g)
        self.assertEqual(a['Flora alba'], 'both')
        self.assertEqual(a['Flora rubra'], 'butterfly_recorded_only')
        self.assertNotIn('Flora missing', a)

    def test_coverage_no_depth_imputation(self):
        df = pd.DataFrame([
            ('Flora alba', 'flower_primary_color', 'white'),
            ('Flora alba', 'floral_form', 'tubular'),
            ('Flora rubra', 'flower_primary_color', 'red_pink'),
            ('Flora rubra', 'floral_form', 'funnel_trumpet'),
            ('Flora rubra', 'tube_depth_class', 'deep'),
            ('Flora longa', 'flower_primary_color', 'white'),
            ('Flora longa', 'floral_form', 'salverform'),
            ('Flora longa', 'tube_depth_class', 'deep'),
            ('Flora mixa', 'flower_primary_color', 'white|red_pink'),
            ('Flora mixa', 'floral_form', 'tubular'),
            ('Flora mixa', 'tube_depth_class', 'deep'),
        ], columns=['accepted_species', 'trait_name', 'normalized_value'])
        status = pd.DataFrame([
            ('A', 'Flora alba', 'introduced'),
            ('B', 'Flora longa', 'introduced'),
            ('B', 'Flora rubra', 'native'),
            ('B', 'Flora mixa', 'unresolved')
        ], columns=['island_id', 'accepted_species', 'origin_status'])
        out = mod.report_scope(df, {'Flora longa': 'moth_recorded_only'}, status, {'A', 'B'})
        self.assertEqual(out['white_single_form_candidate_names'], 2)
        self.assertEqual(out['white_single_form_depth_unknown'], 1)
        self.assertEqual(out['white_single_form_depth_deep'], 1)
        self.assertEqual(out['white_form_guild_positive_evidence']['moth_recorded_only'], 1)
        self.assertEqual(out['white_form_flora']['introduced_pairs_NOT_naturalization'], 2)
        self.assertEqual(out['strict_deep_flora']['islands_with_introduced_record'], 1)

    def test_duplicate_species_trait_rejected(self):
        a = pd.DataFrame([
            ('A a', 'floral_form', 'tubular'), ('A a', 'floral_form', 'salverform')
        ], columns=['accepted_species', 'trait_name', 'normalized_value'])
        with self.assertRaisesRegex(ValueError, 'Duplicated'):
            mod.read_species_traits(a)


if __name__ == '__main__':
    unittest.main()
