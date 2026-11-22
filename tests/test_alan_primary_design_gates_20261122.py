"""Fail-closed validation of source-driven 28-report ALAN triage.

These tests validate provenance and accounting, NOT any new estimated effect size.
"""
import csv
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / 'database' / 'alan_experiments' / 'arthropod_expansion_20261122'


def read_rows(name):
    with (HERE / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


class OriginalSourceScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = read_rows('fulltext_candidates.csv')
        cls.audited = read_rows('primary_design_gates_20261122.csv')
        cls.d = {x['record_id']: x for x in cls.audited}
        cls.summary = json.loads((HERE / 'primary_design_gates_summary_20261122.json').read_text('utf-8'))

    def test_all_28_records_and_order_preserved(self):
        self.assertEqual(len(self.original), 28)
        self.assertEqual([r['record_id'] for r in self.original], [r['record_id'] for r in self.audited])
        self.assertEqual(len(self.d), 28)

    def test_no_new_unearned_meta_effects_or_human_screening(self):
        self.assertTrue(all(r['new_meta_effect_adopted'] == 'false' for r in self.audited))
        self.assertTrue(all(r['raw_data_retrieved_in_this_pass'] == 'false' for r in self.audited))
        self.assertEqual(self.summary['new_eligible_effects'], 0)
        self.assertIsNone(self.summary['fully_independent_eligible_programmes'])
        self.assertFalse(self.summary['human_double_screening'])

    def test_verified_identity_corrections_only(self):
        known = {
            'AGR:IND609454928':'10.1016/j.biocon.2026.111711',
            'AGR:IND609395919':'10.1016/j.biocon.2025.111638',
            'AGR:IND609426870':'10.1016/j.biocon.2026.111722',
            'AGR:IND609249779':'10.1016/j.biocon.2025.111390',
        }
        for k,v in known.items():
            self.assertEqual(self.d[k]['canonical_doi'],v)
            self.assertEqual(self.d[k]['original_doi'],'')

    def test_same_cricket_preprint_not_new_k(self):
        final = self.d['MED:38531451']
        preprint = self.d['PPR:PPR753662']
        self.assertEqual(final['programme_id'], preprint['programme_id'])
        self.assertEqual(preprint['design_gate'],'duplicate_preprint_final')
        self.assertEqual(len({r['programme_id'] for r in self.audited}),27)

    def test_other_programme_overlap_and_stops(self):
        swiss = self.d['MED:38795821']
        lizard = self.d['MED:28139040']
        self.assertEqual(swiss['programme_id'], 'Swiss_ALANex_36lit6dark_3forest')
        self.assertEqual(lizard['design_gate'], 'existing_shared_programme_not_new_k')
        self.assertIn('32 plates', self.d['MED:31284209']['unresolved_gate'])
        self.assertIn('10.34894/SNKXWK', self.d['AGR:IND609454928']['unresolved_gate'])

    def test_counts_match_and_no_false_primary_effect(self):
        from collections import Counter
        c = dict(Counter(x['design_gate'] for x in self.audited))
        self.assertEqual(c, self.summary['counts_by_gate'])
        self.assertEqual(c['core_capture_candidate_variance_hold'],3)
        self.assertEqual(c['separate_attraction_quantitative_candidate'],4)
        self.assertEqual(c['existing_shared_programme_not_new_k'],1)
        self.assertTrue(all(x['evidence_level'] and x['primary_source_url'] and x['unresolved_gate'] for x in self.audited))


if __name__ == '__main__':
    unittest.main()
