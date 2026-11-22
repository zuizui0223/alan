"""Fail-closed original-source and record-grain gates for ALAN expansion.

No review counts are interpreted as independent experiments; this checks the
2026-11-22 search snapshot and the exact 65-item imported source archive.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'database/alan_experiments/arthropod_expansion_20261122'


def csv_rows(relative_path: str) -> list[dict[str, str]]:
    with (BASE / relative_path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def load_json(relative_path: str):
    return json.loads((BASE / relative_path).read_text(encoding='utf-8'))


class TestAlanArthropod20261122(unittest.TestCase):
    def test_exact_source_file_sha256(self):
        checks = load_json('package_files_sha256.json')
        self.assertEqual(len(checks), 65)
        for relative, digest in checks.items():
            with self.subTest(path=relative):
                path = ROOT / relative
                self.assertTrue(path.is_file())
                self.assertTrue(path.resolve().is_relative_to(ROOT.resolve()))
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)

    def test_screening_distinguishes_records_from_experiments(self):
        records = csv_rows('records.csv')
        screened = csv_rows('screened_unmatched_records.csv')
        fulltexts = csv_rows('fulltext_candidates.csv')
        identity = csv_rows('identity_candidates.csv')
        self.assertEqual((len(records), len(screened), len(fulltexts), len(identity)), (341, 143, 28, 15))
        self.assertEqual(len({x['record_id'] for x in records}), 341)
        self.assertEqual(len({x['record_id'] for x in screened}), 143)
        self.assertEqual({x['record_id'] for x in fulltexts},
                         {x['record_id'] for x in screened if x['decision'] == 'fulltext_needed'})
        statuses = Counter(x['decision'] for x in screened)
        self.assertEqual(statuses, {
            'fulltext_needed': 28,
            'identity_review': 7,
            'out_of_scope': 56,
            'reference_mining': 20,
            'separate_indoor': 9,
            'separate_observational': 21,
            'separate_spectrum_no_dark': 2,
        })
        summary = load_json('screening_summary.json')
        self.assertEqual(summary['screening_decisions'], statuses)
        self.assertFalse(summary['final_eligibility'])
        self.assertIsNone(summary['independent_experiment_count'])
        self.assertEqual(summary['identity_candidate_pairs'], 15)

    def test_ris_and_reference_trail_not_inflated(self):
        lines = (BASE / 'unmatched_records.ris').read_text(encoding='utf-8').splitlines()
        self.assertEqual(sum(line.startswith('TY  -') for line in lines), 143)
        self.assertEqual(sum(line.startswith('ER  -') for line in lines), 143)
        citations = csv_rows('citation_ledger.csv')
        self.assertEqual(len(citations), 49)
        self.assertEqual(len({row['report_key'] for row in citations}), 49)
        outcomes = Counter(row['incremental_status'] for row in citations)
        self.assertEqual(outcomes, {'additional_citation_candidate': 31, 'in_expansion_search_only': 9, 'in_prior_inventory': 9})
        manifest = load_json('search_manifest.json')
        self.assertEqual(manifest['raw_hits'], 341)
        self.assertEqual(manifest['records_matching_baseline'], 198)
        self.assertFalse(manifest['final_prisma'])

    def test_original_capture_covariance_and_independence(self):
        estimates = {row['arm']: row for row in csv_rows('data_sources/robert2025/paired_capture_checks.csv')}
        cov = {(row['treatment_a'],row['treatment_b']):float(row['sample_covariance_of_paired_estimates'])
               for row in csv_rows('data_sources/robert2025/paired_capture_estimate_covariance.csv')}
        arms = ('A', 'FW', 'CW')
        self.assertEqual(len(cov),9)
        for i in arms:
            self.assertEqual(int(estimates[i]['n_paired_reserves']),4)
            self.assertTrue(math.isclose(cov[i,i],float(estimates[i]['paired_difference_in_changes_se'])**2,
                                         rel_tol=1e-10,abs_tol=1e-12))
            for j in arms:
                self.assertTrue(math.isclose(cov[i,j],cov[j,i],rel_tol=1e-12,abs_tol=1e-12))
        self.assertAlmostEqual(float(estimates['A']['paired_difference_in_changes_mean']),3.583333333333333,places=9)
        self.assertAlmostEqual(float(estimates['FW']['paired_difference_in_changes_mean']),9.083333333333334,places=9)
        self.assertAlmostEqual(float(estimates['CW']['paired_difference_in_changes_mean']),14.083333333333332,places=9)
        # Do not treat the all-lit alternative summary as an independent fourth effect.
        self.assertEqual(len(estimates),4)
        self.assertIn('all_lit_site_mean',estimates)

    def test_swiss_source_overlap_and_flight_stop(self):
        overlaps = load_json('data_sources/swiss_envidat/ground_family_overlap.json')
        self.assertEqual(overlaps['eligible_566_ground_rows'],307)
        self.assertEqual(overlaps['exact_matches'],306)
        self.assertEqual(len(overlaps['unmatched']),1)
        flagged = csv_rows('data_sources/swiss_envidat/flagged_flight_records.csv')
        self.assertGreaterEqual(len(flagged),28)
        source = load_json('data_sources/swiss_envidat/source_manifest.json')
        self.assertEqual(len(source['resources']),14)
        for row in source['resources']:
            file = BASE / 'data_sources/swiss_envidat' / row['dataset'] / row['filename']
            self.assertTrue(file.is_file())
            self.assertEqual(file.stat().st_size,row['bytes'])
            self.assertEqual(hashlib.sha256(file.read_bytes()).hexdigest(),row['sha256'])

    def test_no_false_external_search_or_pooled_result(self):
        access = csv_rows('contract_database_access_20261122.csv')
        self.assertEqual({x.get('database') for x in access},{'Web of Science Core Collection','Scopus'})
        checks = load_json('verification_summary.json')
        self.assertFalse(checks['new_pooled_estimate'])
        self.assertFalse(checks['human_double_extraction_completed'])
        self.assertEqual(checks['existing_effect_catalogue_rows'],108)
        self.assertEqual(checks['contract_database_browser_attempts']['queries_submitted'],0)


if __name__ == '__main__':
    unittest.main()