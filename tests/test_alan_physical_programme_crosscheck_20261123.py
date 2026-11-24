import csv, collections, unittest
from pathlib import Path
PATH = Path(__file__).resolve().parents[1] / 'database/alan_experiments/arthropod_expansion_20261122/cross_programme/physical_programme_ledger.csv'
class ProgrammeIdentityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  with PATH.open(encoding='utf-8', newline='') as f:cls.rows=list(csv.DictReader(f))
 def test_seven_distinct_physical_programmes(self):
  # Exact physical programme identities; two or three papers per programme never
  # multiply independent lamps, flumes, or the metanalytic study denominator.
  expected={
   'Michigan2016_Parkinson20','Ohio2013_Sullivan15',
   'Westhavelland2012_2014_Manfrin',
   'Landau_RSM_2022_16flumes_ALAN_crayfish',
   'SwissGurbe2021_Carannante',
   'EERES2021_Disconnected18flumes_6headers',
   'SevenPonds2022_Parkinson30plots'
  }
  self.assertEqual({r['programme_id'] for r in self.rows},expected)
 def test_eeres_two_publications_share_one_original_six_header_programme(self):
  e=[r for r in self.rows if r['programme_id']=='EERES2021_Disconnected18flumes_6headers']
  self.assertEqual({r['doi'] for r in e},{'10.1002/ecog.07758','10.1016/j.scitotenv.2023.168836'})
  self.assertEqual({r['experiment_year'] for r in e},{'2021'})
  self.assertTrue(all(r['meta_status']!='newly_approved_meta_effect' for r in e))
  self.assertTrue(all('6_header_tanks_3_flumes_each' in r['physical_light_allocation'] for r in e))
 def test_sevenponds_2023_and_2026_share_2022_lamps(self):
  p=[r for r in self.rows if r['programme_id']=='SevenPonds2022_Parkinson30plots']
  self.assertEqual({r['doi'] for r in p},{'10.1098/rstb.2022.0361','10.1111/fwb.70236'})
  self.assertEqual({r['experiment_year'] for r in p},{'2022'})
  self.assertTrue(all('30_pond_margin' in r['physical_light_allocation'] for r in p))
 def test_eeres_and_landau_are_distinct_despite_shared_authors(self):
  e=[r for r in self.rows if r['programme_id'].startswith('EERES2021_')]
  l=[r for r in self.rows if r['programme_id'].startswith('Landau_RSM_2022_')]
  self.assertEqual({r['experiment_year'] for r in e},{'2021'})
  self.assertEqual({r['experiment_year'] for r in l},{'2022'})
 def test_landau_three_papers_one_programme(self):
  p=[r for r in self.rows if r['programme_id']=='Landau_RSM_2022_16flumes_ALAN_crayfish']
  self.assertEqual(len({r['doi'] for r in p}), 3)
  self.assertEqual({r['experiment_year'] for r in p}, {'2022'})
 def test_westhavelland_distinct_from_landau(self):
  w=[r for r in self.rows if r['programme_id'].startswith('Westhavelland')]
  l=[r for r in self.rows if r['programme_id'].startswith('Landau')]
  self.assertNotEqual({r['programme_id'] for r in w},{r['programme_id'] for r in l})
 def test_gurbe_original_field_year_2017(self):
  rows=[r for r in self.rows if r['programme_id']=='SwissGurbe2021_Carannante']
  self.assertEqual({r['experiment_year'] for r in rows},{'2017'})
 def test_no_new_approved_meta(self):
  self.assertTrue(all(r['meta_status']!='newly_approved_meta_effect' for r in self.rows))
 def test_response_types_not_pooled(self):
  m=[r for r in self.rows if r['programme_id']=='Michigan2016_Parkinson20']
  self.assertEqual(len({r['measurement'] for r in m}),3)
 def test_undereplicated_programme_hold(self):
  w=[r for r in self.rows if r['programme_id'].startswith('Westhavelland')]
  self.assertTrue(all(r['meta_status']=='descriptive_only' for r in w))
if __name__=='__main__':unittest.main()
