"""Originally fetched official-source manifest checks only, no biological replications."""
import json
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'database/alan_experiments/arthropod_expansion_20261122/public_source_receipts'
class PublicSourceReceiptsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.j=json.loads((ROOT/'charvalakis_reinhard_original_source_gate_20261122.json').read_text())
        cls.b=json.loads((ROOT/'carannante_boris_access_gate_20261122.json').read_text())
    def test_dataverse_is_restricted_not_missing_dataset(self):
        s=self.j['sources']['charvalakis2026'];self.assertTrue(s['meta_verified'])
        self.assertEqual(s['file_count'],4);self.assertEqual(s['downloaded_tabular_files'],0)
        self.assertEqual(sum(bool(x['restricted']) for x in s['file_registry']),4)
    def test_roadlight_raw_counts_and_real_fixture_n(self):
        reg=self.j['sources']['reinhard2026']['file_registry']
        x=next(x for x in reg if x['label']=='Roadlights_long.csv')
        self.assertEqual((x['row_count'],x['column_count'],x['delimiter']),(701,26,';'))
        self.assertEqual(x['download_status'],'ORIGINAL_BYTES_VERIFIED')
        g=x['raw_count_by_treatment_method'];self.assertEqual(len(g),2)
        self.assertEqual([v['raw_count_column_sum'] for v in g],[1198.,511.])
        self.assertEqual([v['n_trap_labels'] for v in g],[3,3]);self.assertEqual([v['n_date_labels'] for v in g],[13,13])
        self.assertEqual([v['n_interval_labels'] for v in g],[5,5])
    def test_source_is_no_dark_not_new_meta_k(self):
        s=self.j['sources']['reinhard2026']
        self.assertIn('no unlit control',s['endpoint'])
        self.assertEqual(s['effect_sizes_estimated'],0)
        self.assertEqual(s['independent_study_k_increase'],0)
    def test_boris_redirect_to_original_dataset_item_but_no_raw_rows(self):
        b=self.b
        self.assertFalse(b['original_plot_by_night_numeric_effects_extracted'])
        self.assertEqual(b['independent_programmes_added'],0)
        self.assertEqual(b['probes'][0]['final_url'],'https://boris-portal.unibe.ch/items/f559ab1e-e693-4fbc-b7fe-683134f6dbe1')
        self.assertTrue(all(x['status']=='PUBLIC_HTML_RECEIVED' for x in b['probes']))
if __name__=='__main__':unittest.main()
