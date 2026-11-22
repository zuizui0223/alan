"""Independent arithmetic audit of frozen 6-header EERES aggregates and joint covariance."""
import json,math,statistics,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'database/alan_experiments/arthropod_expansion_20261122/eeres2021_source_reanalysis/original_reanalysis_20261123.json'
sys.path.insert(0,str(ROOT/'scripts/analysis'))
from alan_eeres2021_verify_frozen_20261123 import verify

class FrozenEERES(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.a=json.loads(SRC.read_text(encoding='utf-8'))
 def test_source_sha_and_unit(self):
  self.assertEqual(self.a['source_sha256'],'496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f')
  self.assertEqual(self.a['joint_shared_HT']['independent_header_blocks'],6)
  self.assertEqual(len(self.a['joint_shared_HT']['source_outcomes']),8)
 def test_original_8_effects_independent_replay(self):
  for endpoint,d in self.a['effect_summaries'].items():
   for metric,entries in d.items():
    self.assertEqual(len(entries),2)
    for e in entries:
     self.assertEqual(len(e['block_treatment_values']),6)
     aa=statistics.mean(e['block_treatment_values'])
     bb=statistics.mean(e['block_control_values'])
     self.assertAlmostEqual(e['lnRR'],math.log(aa/bb),places=12)
     self.assertAlmostEqual(e['RR'],aa/bb,places=12)
     vari=statistics.variance(e['block_treatment_values'])
     varj=statistics.variance(e['block_control_values'])
     cov=statistics.covariance(e['block_treatment_values'],e['block_control_values'])
     var=(vari/aa**2+varj/bb**2-2*cov/aa/bb)/6
     self.assertAlmostEqual(e['SE_paired_HT_delta'],math.sqrt(max(0,var)),places=12)
 def test_correlation_aware_covariance_and_contrasts(self):
  r=self.a['joint_shared_HT'];c=r['covariance_logRR'];labels=r['source_outcomes']
  self.assertEqual(len(c),8)
  estimates={}
  for endpoint,d in self.a['effect_summaries'].items():
   for metric,entries in d.items():
    for e in entries:
     key=f'{endpoint}|{metric}|{"L_vs_CT" if e["treatment"]=="L" else "LF_vs_F"}'
     estimates[key]=e
  for i,ki in enumerate(labels):
   self.assertEqual(len(c[i]),8)
   self.assertAlmostEqual(c[i][i],estimates[ki]['SE_paired_HT_delta']**2,places=12)
   for j,kj in enumerate(labels):self.assertAlmostEqual(c[i][j],c[j][i],places=12)
  for q in r['stage_contrasts']+r['light_by_flow_modifiers']:
   i=labels.index(q['numerator_lnRR']);j=labels.index(q['denominator_lnRR'])
   expected=math.exp(estimates[labels[i]]['lnRR']-estimates[labels[j]]['lnRR'])
   se=math.sqrt(max(0,c[i][i]+c[j][j]-2*c[i][j]))
   self.assertAlmostEqual(q['ratio_of_response_ratios'],expected,places=12)
   self.assertAlmostEqual(q['SE_shared_HT'],se,places=12)
 def test_important_stage_contrasts(self):
  c={r['numerator_lnRR']:r for r in self.a['joint_shared_HT']['stage_contrasts']}
  self.assertAlmostEqual(c['drift|abund|L_vs_CT']['ratio_of_response_ratios'],1.668980606126087,places=12)
  self.assertAlmostEqual(c['drift|abund|LF_vs_F']['ratio_of_response_ratios'],0.7004766369745273,places=12)
  mods={r['numerator_lnRR']:r for r in self.a['joint_shared_HT']['light_by_flow_modifiers']}
  self.assertAlmostEqual(mods['benthos|abund|LF_vs_F']['ratio_of_response_ratios'],1.6761141374297504,places=12)
 def test_source_observed_stop_not_fake_zero(self):
  observed={'status':'SOURCE_STOP_ARCHIVE','original_workbook':{'sha256':self.a['source_sha256']}}
  with self.assertRaises((AssertionError,KeyError)):verify(observed,self.a)

if __name__=='__main__':unittest.main()
