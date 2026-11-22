import csv, collections, unittest
from pathlib import Path
PATH = Path(__file__).resolve().parents[1] / 'database/alan_experiments/arthropod_expansion_20261122/cross_programme/physical_programme_ledger.csv'
class ProgrammeIdentityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  with PATH.open(encoding='utf-8', newline='') as f:cls.rows=list(csv.DictReader(f))
 def test_five_real_programmes(self):
  self.assertEqual(len({r['programme_id'] for r in self.rows}), 5)
 def test_landau_three_papers_one_programme(self):
  p=[r for r in self.rows if r['programme_id']=='Landau_RSM_2022_16flumes_ALAN_crayfish']
  self.assertEqual(len({r['doi'] for r in p}), 3)
  self.assertEqual({r['experiment_year'] for r in p}, {'2022'})
 def test_westhavelland_distinct_from_landau(self):
  w=[r for r in self.rows if r['programme_id'].startswith('Westhavelland')]
  l=[r for r in self.rows if r['programme_id'].startswith('Landau')]
  self.assertNotEqual({r['programme_id'] for r in w},{r['programme_id'] for r in l})
 def test_no_new_approved_meta(self):
  self.assertTrue(all(r['meta_status']!='newly_approved_meta_effect' for r in self.rows))
 def test_response_types_not_pooled(self):
  m=[r for r in self.rows if r['programme_id']=='Michigan2016_Parkinson20']
  self.assertEqual(len({r['measurement'] for r in m}),3)
 def test_undereplicated_programme_hold(self):
  w=[r for r in self.rows if r['programme_id'].startswith('Westhavelland')]
  self.assertTrue(all(r['meta_status']=='descriptive_only' for r in w))
if __name__=='__main__':unittest.main()
