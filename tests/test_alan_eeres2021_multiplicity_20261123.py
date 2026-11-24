import json,math,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/analysis'))
from alan_eeres2021_multiplicity_audit_20261123 import t5_two_sided_p,corrections,derive
FILE=ROOT/'database/alan_experiments/arthropod_expansion_20261122/eeres2021_source_reanalysis/original_reanalysis_20261123.json'
class PosthocChecks(unittest.TestCase):
 def test_analytic_t5_known(self):
  self.assertAlmostEqual(t5_two_sided_p(0.,1.),1.,places=12)
  self.assertAlmostEqual(t5_two_sided_p(2.,1.),0.10193947882985836,places=12)
  self.assertAlmostEqual(t5_two_sided_p(10.,1.),0.00017094757574295794,places=12)
 def test_invalid_se_stops(self):
  with self.assertRaises(ValueError):t5_two_sided_p(1,0)
 def test_bh_holm_monotone(self):
  p=[0.01,0.04,0.03,1.0];h,b=corrections(p)
  self.assertEqual(h,[0.04,0.09,0.09,1.0])
  self.assertAlmostEqual(b[0],0.04)
 def test_8_plus_8_posthoc(self):
  x=derive(json.loads(FILE.read_text(encoding='utf8')))
  self.assertEqual(len(x),16)
  self.assertEqual(sum(y['family']=='original_8_effects' for y in x),8)
  self.assertEqual(sum(y['family']=='joint_8_posthoc_contrasts' for y in x),8)
  # Only one original family effect survives 8-test Holm; no joint contrast survives BH.
  self.assertEqual(sum(y['p_holm_within_8']<0.05 for y in x[:8]),1)
  self.assertEqual(sum(y['q_BH_within_8']<0.05 for y in x[8:]),0)
  self.assertAlmostEqual(min(y['q_BH_within_8'] for y in x[8:]),0.081135,delta=1e-5)
if __name__=='__main__':unittest.main()
