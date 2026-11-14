import importlib.util
import math
from pathlib import Path
import unittest
P=Path(__file__).resolve().parents[1]/'scripts/analysis/alan_reproductive_meta_20261114.py'
s=importlib.util.spec_from_file_location('meta',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def row(study,g,se,endpoint='seed set'):
    return {'Study':study,'DOI':study,'SMD':str(g),'SE':str(se),'Measurement':endpoint}
class ReproductiveMeta(unittest.TestCase):
    def test_correlated_mean_variance(self):
        a=m.study_effects([row('a',1,1),row('a',3,1)],.5)
        self.assertEqual(len(a),1)
        self.assertAlmostEqual(a[0]['yi'],2)
        self.assertAlmostEqual(a[0]['vi'],.75)
    def test_perfect_correlation_does_not_manufacture_precision(self):
        a=m.study_effects([row('a',1,.3)]*100,1)
        self.assertAlmostEqual(a[0]['vi'],.09)
        self.assertEqual(a[0]['n_source_effects'],100)
    def test_distinct_endpoints_rejected(self):
        with self.assertRaises(ValueError):m.study_effects([row('a',1,.3),row('b',1,.3,'fruit set')],.5)
    def test_bad_standard_error_rejected(self):
        for se in [0,-1,float('nan')]:
            with self.assertRaises(ValueError):m.study_effects([row('a',1,se)],.5)
    def test_too_few_studies_not_promoted_to_global_result(self):
        z=m.random_effects([{'yi':1,'vi':.1},{'yi':2,'vi':.2}])
        self.assertEqual(z['status'],'TOO_FEW_STUDIES_FOR_POOLED_INFERENCE')
    def test_identical_effects_have_zero_heterogeneity(self):
        z=m.random_effects([{'yi':1,'vi':.1}]*4)
        self.assertAlmostEqual(z['pooled_SMD'],1)
        self.assertAlmostEqual(z['tau2'],0,places=6)
        self.assertGreater(z['ci95'][1],1)
        self.assertLess(z['ci95'][0],1)
if __name__=='__main__':unittest.main()
