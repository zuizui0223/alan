import sys, unittest, math, itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'/'analysis'))
import alan_landau16_original_rds_20261123 as a

class LandauOriginalSourceTests(unittest.TestCase):
 def mock(self):
  d=[]
  for k,t in enumerate(a.TREATMENTS):
   for f in range(4):
    for w,v in [('Week 1',1.0),('Week 4',2.0)]:
     d.append({'treat':t,'week':w,'flume':str(k*4+f+1),'abund_rate':str((k+1)*v)})
  return d
 def test_pinned_original_paths(self):
  self.assertEqual(len(a.SOURCES),2)
  self.assertTrue(a.source_url(a.SOURCES['emergence']['path']).startswith('https://raw.githubusercontent.com/colinog/'))
  with self.assertRaises(ValueError):a.source_url('some_other_source.rds')
 def test_git_sha1_bytes(self):
  self.assertEqual(a.git_blob_sha(b'hello'),'b6fc4c620b67d95f953a5c1c1230aaab5db5a1b0')
 def test_all_physical_flumes_present(self):
  q=a.primary_results(self.mock())
  self.assertEqual(q['grain']['n_real_flumes'],16)
  self.assertEqual(len(q['contrasts']),10)
  self.assertEqual(q['grain']['treatments'],dict.fromkeys(a.TREATMENTS,4))
 def test_physical_flume_duplicate_rejected(self):
  dat=self.mock();dat.append(dict(dat[0]))
  with self.assertRaisesRegex(ValueError,'duplicate'):a.primary_results(dat)
 def test_wrong_unit_label_rejected(self):
  dat=self.mock();dat[0]['flume']='100'
  with self.assertRaisesRegex(ValueError,'16 original flumes'):a.primary_results(dat)
 def test_paired_weeks_not_independent_k(self):
  r=a.primary_results(self.mock());self.assertEqual(r['grain']['weeks'],{'Week 1':16,'Week 4':16})
  self.assertEqual(r['contrasts'][0]['n_treatment_flumes'],4)
 def test_single_epoch_arithmetic_rr(self):
  r=a.primary_results(self.mock());first=r['contrasts'][0]
  self.assertEqual(first['contrast'],'ALAN_no_cray');self.assertAlmostEqual(first['RR'],2)
  self.assertEqual(first['permutation_allocations'],70)
 def test_factorial_interaction(self):
  r=a.primary_results(self.mock());first=[x for x in r['contrasts'] if x['contrast']=='factorial_ALAN_x_Cray_interaction'][0]
  self.assertAlmostEqual(first['RR'],(4/3)/(2/1))
 def test_zero_arm_mean_stop(self):
  d=self.mock()
  for r in d:
   if r['treat']=='control':r['abund_rate']='0'
  with self.assertRaisesRegex(ValueError,'zero arm mean'):a.primary_results(d)
 def test_source_alias_is_one_programme(self):
  self.assertEqual(a.SOURCE_COMMIT,'5e749c03161ad9ae7d57f12ddb86d9aaad283077')

if __name__=='__main__':unittest.main()
