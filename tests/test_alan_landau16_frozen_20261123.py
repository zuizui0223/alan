"""Independent arithmetic replay against source-derived 32 true flume-week means.
The original author RDS is SHA-checked by the separate source CI; no raw records
are added to ALAN. This test verifies frozen derived effects and unit counts.
"""
import collections, csv, itertools, json, math, statistics, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'database/alan_experiments/arthropod_expansion_20261122/landau2022_source'
JSON_FILE=BASE/'original_reanalysis_20261123.json'
CSV_FILE=BASE/'effect_candidates_20261123.csv'

def ratio_and_perm(x,y):
  x=list(x);y=list(y)
  actual=math.log(statistics.mean(x)/statistics.mean(y))
  samples=x+y;count=0;valid=0
  for ind in itertools.combinations(range(8),4):
    ids=set(ind)
    a=[v for i,v in enumerate(samples) if i in ids]
    b=[v for i,v in enumerate(samples) if i not in ids]
    if min(statistics.mean(a),statistics.mean(b))<=0:continue
    curr=abs(math.log(statistics.mean(a)/statistics.mean(b)))
    valid+=1
    if curr+1e-12>=abs(actual):count+=1
  return math.exp(actual),count/valid,valid

class LandauFrozenSourceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.doc=json.loads(JSON_FILE.read_text(encoding='utf-8'))
  with CSV_FILE.open(newline='',encoding='utf-8') as f: cls.effects=list(csv.DictReader(f))
 def test_pinned_original_data_and_one_programme(self):
  x=self.doc
  self.assertEqual(x['programme_id'],'Landau_RSM_2022_16flumes_ALAN_crayfish')
  sources={z['path']:z for z in x['original_source_receipts']}
  self.assertEqual(sources['data/other data_output/emergence.rds']['sha256'],'289e9dc874724ca765aa2e8314f975ef7cf7b2cb688bdb62e8f1be9bc50d1776')
  self.assertEqual(sources['data/other data_output/other_emergence.rds']['sha256'],'6723d8435bab80dca1d73d8029b5a32f16dff4689923bcde4185dfb01b4d8dc0')
  self.assertEqual(x['number_new_independent_programmes'],0)
 def test_source_physical_unit(self):
  g=self.doc['original_primary']['grain']
  self.assertEqual(g['n_rows'],192)
  self.assertEqual(g['n_real_flumes'],16)
  self.assertEqual(g['n_flume_week_means'],32)
  self.assertEqual(g['subsamples_per_flume_week'],6)
  self.assertEqual(sorted(g['treatments'].values()),[4,4,4,4])
  self.assertEqual(self.doc['other_emergence_grain']['n_rows'],529)
  self.assertTrue(self.doc['other_emergence_grain']['unbalanced_taxon_rows_preclude_straight_sum'])
 def test_no_flume_pseudoreplication(self):
  rows=self.doc['original_primary']['aggregated_flume_week_means']
  self.assertEqual(len(rows),32)
  self.assertEqual(len({r['flume'] for r in rows}),16)
  self.assertEqual({r['week'] for r in rows},{'Week 1','Week 4'})
  self.assertTrue(all(r['n_trap_day_subsamples']==6 for r in rows))
  self.assertEqual(len({(r['flume'],r['week']) for r in rows}),32)
  assignments=collections.defaultdict(set)
  for r in rows:assignments[r['flume']].add(r['treatment'])
  self.assertTrue(all(len(v)==1 for v in assignments.values()))
 def test_all_8_arm_vs_arm_effects_from_flume_values(self):
  obs=self.doc['original_primary']['aggregated_flume_week_means']
  got=0
  for r in self.doc['original_primary']['contrasts']:
   if r['contrast']=='factorial_ALAN_x_Cray_interaction':continue
   cells={arm:[z['mean_abund_rate'] for z in obs if z['week']==r['week'] and z['treatment']==arm] for arm in ['alan','control','cray','alan+cray']}
   a,b=cells[r['treatment']],cells[r['reference']]
   self.assertEqual((len(a),len(b)),(4,4))
   rr,p,n=ratio_and_perm(a,b)
   self.assertAlmostEqual(rr,r['RR'],places=10)
   self.assertAlmostEqual(p,r['exact_4of8_label_permutation_two_sided_p'],places=12)
   self.assertEqual(n,70)
   v=statistics.variance(a)/(4*statistics.mean(a)**2)+statistics.variance(b)/(4*statistics.mean(b)**2)
   self.assertAlmostEqual(math.sqrt(v),r['delta_SE'],places=10)
   got+=1
  self.assertEqual(got,8)
 def test_factorial_interactions_reconstructed(self):
  for week in ('Week 1','Week 4'):
   effects={r['contrast']:r for r in self.doc['original_primary']['contrasts'] if r['week']==week}
   expected=effects['ALAN_with_cray']['RR']/effects['ALAN_no_cray']['RR']
   self.assertAlmostEqual(effects['factorial_ALAN_x_Cray_interaction']['RR'],expected,places=10)
 def test_two_repeated_week_contrasts_not_new_k(self):
  source=self.doc['original_primary']['aggregated_flume_week_means']
  by_flume=collections.defaultdict(dict)
  for r in source:by_flume[r['flume']][r['week']]=(r['treatment'],r['mean_abund_rate'])
  for r in self.doc['original_primary']['temporal_within_flume_change']:
   light,reference=('alan','control') if r['contrast']=='ALAN_no_cray' else ('alan+cray','cray')
   vals={light:[],reference:[]}
   for f,wks in by_flume.items():
    trt,x=wks['Week 1'];_,y=wks['Week 4']
    if trt in vals:vals[trt].append(math.log(y/x))
   self.assertEqual((len(vals[light]),len(vals[reference])),(4,4))
   delta=statistics.mean(vals[light])-statistics.mean(vals[reference])
   self.assertAlmostEqual(math.exp(delta),r['ratio_of_geometric_changes'],places=10)
   obs=abs(delta);allvals=vals[light]+vals[reference];nn=0
   for ids in itertools.combinations(range(8),4):
    aa=set(ids);a=[v for k,v in enumerate(allvals) if k in aa];b=[v for k,v in enumerate(allvals) if k not in aa]
    if abs(statistics.mean(a)-statistics.mean(b))+1e-12>=obs:nn+=1
   self.assertAlmostEqual(nn/70,r['exact_4of8_two_sided_label_permutation_p'],places=12)
 def test_effect_csv_matches_source_output(self):
  self.assertEqual(len(self.effects),10)
  self.assertEqual({r['effect_status'] for r in self.effects},{'exploratory_not_canonical'})
  for row in self.effects:
   got=[r for r in self.doc['original_primary']['contrasts'] if r['week']==row['week'] and r['contrast']==row['contrast']]
   self.assertEqual(len(got),1)
   self.assertAlmostEqual(float(row['RR']),got[0]['RR'],places=9)

if __name__=='__main__':unittest.main()
