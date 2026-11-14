import importlib.util
import math
from pathlib import Path
import unittest
P=Path(__file__).resolve().parents[1]/'scripts/analysis/alan_spillover_effects_20261114.py'
spec=importlib.util.spec_from_file_location('audit',P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class Effects(unittest.TestCase):
    def test_log_odds_interval(self):
        r=m.effect(1.90, .576)
        self.assertAlmostEqual(r['ratio'], math.exp(1.90))
        self.assertAlmostEqual(r['log_variance'], .576**2)
        self.assertLess(r['ratio_ci95'][0],r['ratio'])
    def test_reverse_contrast(self):
        a=m.effect(1.90,.576); b=m.effect(-1.90,.576)
        self.assertAlmostEqual(a['ratio']*b['ratio'],1)
        self.assertAlmostEqual(a['ratio_ci95'][0]*b['ratio_ci95'][1],1)
    def test_no_invalid_standard_error(self):
        for se in [0,-1,float('nan')]:
            with self.assertRaises(ValueError):m.effect(1,se)
    def test_table_preserves_cell_labels(self):
        import xml.etree.ElementTree as E
        t=E.fromstring('<table><tr><th>term</th><th>beta</th></tr><tr><td rowspan="2">seed</td><td>1.0</td></tr><tr><td>2.0</td></tr></table>')
        r=m.table_rows(t)
        self.assertEqual(r[1][0]['rowspan'],2)
        self.assertEqual(r[2][0]['text'],'2.0')
    def test_different_outcomes_not_pooled(self):
        with self.assertRaises(ValueError):m.common_endpoint([{'endpoint':'seed'},{'endpoint':'visits'}])
        self.assertEqual(m.common_endpoint([{'endpoint':'seed'},{'endpoint':'seed'}]),'seed')
if __name__=='__main__':unittest.main()
