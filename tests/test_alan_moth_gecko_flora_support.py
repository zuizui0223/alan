import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts' / 'analysis'))
import alan_moth_syndrome_validate as moth
import alan_island_gecko_light_support as gecko


def write_csv(path, columns, rows):
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p = Path(self.tmp.name)
        self.profiles = self.p / 'profiles.csv'
        self.traits = self.p / 'traits.csv'
        self.interactions = self.p / 'globi.csv'
        self.floras = self.p / 'floras.csv'
        self.geckos = self.p / 'geckos.csv'
        self.light = self.p / 'light.csv'
        self.islands = self.p / 'islands.csv'
        write_csv(self.profiles, ['profile_id', 'flower_primary_color', 'floral_form'], [
            [0, 'white', 'tubular'], [1, 'red_pink', 'open_radial'],
            [2, 'white|red_pink', 'tubular'], [3, 'white', 'salverform']])
        write_csv(self.traits, ['accepted_species', 'flower_primary_color',
                                'floral_form', 'tube_depth_class'], [
            ['Flora alba', 'white', 'tubular', 'deep'],
            ['Flora rubra', 'red_pink', 'open_radial', 'shallow'],
            ['Flora mixa', 'white|red_pink', 'tubular', 'deep'],
            ['Flora longa', 'white', 'salverform', 'deep'],
            ['Flora brevis', 'white', 'tubular', 'unresolved'],
        ])
        write_csv(self.interactions, ['plant_name', 'guild', 'relation'], [
            ['Flora alba', 'moth_family_visit', 'visitsFlowersOf'],
            ['Flora rubra', 'butterfly_family_visit', 'visitsFlowersOf'],
            ['Flora longa', 'moth_family_visit', 'pollinates'],
            ['Flora longa', 'butterfly_family_visit', 'visitsFlowersOf'],
            ['Flora mixa', 'moth_family_visit', 'visitsFlowersOf'],
        ])
        write_csv(self.floras, ['island_id', 'accepted_species', 'establishment_status',
                                'review_status', 'establishment_year'], [
            ['A', 'Flora alba', 'naturalized', 'accepted', 2021],
            ['A', 'Flora rubra', 'naturalized', 'accepted', ''],
            ['A', 'Flora mixa', 'native', 'accepted', ''],
            ['B', 'Flora longa', 'naturalized', 'accepted', 2010],
            ['B', 'Flora rubra', 'naturalized', 'accepted', ''],
            ['C', 'Flora alba', 'cultivated', 'accepted', ''],
            ['C', 'Flora longa', 'naturalized', 'accepted', 2020],
        ])
        write_csv(self.geckos, ['island_id', 'gecko_status', 'review_status',
                                'geographic_scope', 'source'], [
            ['A', 'present', 'accepted', 'island', 'reviewed survey'],
            ['B', 'confirmed_absent', 'accepted', 'island', 'reviewed checklist'],
            ['C', 'unknown', 'accepted', 'island', 'occurrence gap'],
        ])
        write_csv(self.light, ['island_id', 'year', 'radiance', 'product_version',
                               'n_valid_months', 'coverage_fraction'], [
            ['A', 2020, 9, 'VNL V2.2', 12, 1], ['A', 2021, 8, 'VNL V2.2', 12, 1],
            ['B', 2020, 1, 'VNL V2.2', 12, 1], ['B', 2021, 2, 'VNL V2.2', 12, 1],
            ['C', 2020, 2, 'VNL V2.2', 12, 1], ['C', 2021, 3, 'VNL V2.2', 12, 1],
        ])
        write_csv(self.islands, ['island_id', 'archipelago', 'area_km2'], [
            ['A', 'Test archipelago', 30], ['B', 'Test archipelago', 40],
            ['C', 'Test archipelago', 80],
        ])

    def test_shape_and_positive_evidence(self):
        r = moth.run(self.profiles, self.traits, self.interactions,
                     self.p / 'shape_result')
        self.assertEqual(r['total_unique_trait_taxa'], 5)
        self.assertEqual(r['counts']['moth_recorded_only']['candidate'], 1)
        self.assertEqual(r['counts']['both']['candidate'], 1)
        self.assertEqual(r['counts']['moth_recorded_only']['unknown'], 1)
        self.assertEqual(r['counts']['no_guild_record']['unknown'], 1)
        self.assertFalse(r['biological_specificity_estimable'])

    def test_island_support(self):
        r = gecko.run(self.profiles, self.traits, self.floras, self.geckos,
                      self.light, self.islands, self.p / 'support', min_naturalized=2)
        self.assertEqual(r['islands_in_reference_metadata'], 3)
        self.assertEqual(r['islands_with_all_three_sources_and_naturalized_threshold'], 2)
        self.assertEqual(r['eligible_archipelagos'], 1)
        self.assertEqual(r['archipelagos_with_both_verified_gecko_statuses'], 1)
        self.assertEqual(r['descriptive_2x2']['present__high']['naturalized_mothlike'], 1)
        self.assertEqual(r['descriptive_2x2']['confirmed_absent__low_or_median']['naturalized_other_known'], 1)
        self.assertEqual(r['eligible_islands_with_dated_post2012_naturalized_records'], 1)

    def test_country_gecko_not_island(self):
        write_csv(self.geckos, ['island_id', 'gecko_status', 'review_status',
                                'geographic_scope', 'source'], [
            ['A', 'present', 'accepted', 'country', 'country list']])
        with self.assertRaisesRegex(ValueError, 'country-only'):
            gecko.load_geckos(self.geckos)

    def test_broad_tier_does_not_impute_tube_depth(self):
        ref = moth.reference_profiles(self.profiles)
        self.assertEqual(moth.broad_phenotype('white', 'tubular', ref), 'candidate')
        self.assertEqual(moth.broad_phenotype('white|red_pink', 'tubular', ref), 'unknown')
        broad = moth.run(self.profiles, self.traits, self.interactions,
                         self.p / 'broad_result', tier='broad')
        self.assertEqual(broad['phenotype_tier'], 'broad')
        self.assertFalse(broad['strict_tier_requires_documented_deep_tube'])
        self.assertEqual(broad['counts']['no_guild_record']['candidate'], 1)

    def test_unknown_not_noncandidate(self):
        ref = moth.reference_profiles(self.profiles)
        self.assertEqual(moth.phenotype('white', 'tubular', 'unresolved', ref), 'unknown')
        self.assertEqual(moth.phenotype('white|red_pink', 'tubular', 'deep', ref), 'unknown')
        self.assertEqual(moth.phenotype('white', 'tubular', 'deep', ref), 'candidate')


if __name__ == '__main__':
    unittest.main()
