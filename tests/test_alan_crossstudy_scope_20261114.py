import math
import unittest
from scripts.analysis.alan_crossstudy_scope_20261114 import parse_references, classify_scope, odds_ratio

class ScopeTests(unittest.TestCase):
    def test_multi_year_one_reference_and_verified_suffix_dedup(self):
        self.assertEqual(parse_references('(Giavi et al., 2021, 2020; Van Langevelde et al., 2017a; Van Langevelde et al., 2017b)'),
                         {'giavi:2020','giavi:2021','van langevelde:2017'})
        self.assertEqual(len(parse_references('(Duque et al., 2021a; Duque et al., 2021b)')),2)
    def test_unicode_author_and_no_record_not_zero(self):
        self.assertEqual(parse_references('(Dzul-Cauich and Munguía-Rosas, 2022; Altermatt and Ebert, 2016)'),
                         {'dzul-cauich:2022','altermatt:2016'})
        with self.assertRaises(ValueError): parse_references('none given')
    def test_reference_membership_is_not_effect_weight(self):
        r=classify_scope(['(A et al., 2000; B et al., 2001)', '(C et al., 2002; A et al., 2000)'],
                         '(B et al., 2001)', '(B et al., 2001; C et al., 2002)')
        self.assertEqual(r['light_papers'],['a:2000','b:2001','c:2002'])
        self.assertEqual(r['alan_candidate_reproductive_papers'],['c:2002'])
        self.assertEqual(r['excluded_uv_reproductive_papers'],['b:2001'])
        self.assertFalse(r['alan_only_effect_reestimated'])
    def test_unknown_exposure_not_included_as_night(self):
        with self.assertRaises(ValueError): classify_scope(['(A et al., 2000)'],'(B et al., 2001)','(A et al., 2000)')
    def test_or_reorientation_no_recruitment_percentage(self):
        r=odds_ratio(-1.04,.499)
        self.assertAlmostEqual(r['odds_ratio'],math.exp(-1.04))
        self.assertLess(r['ci95'][1],1)
        self.assertEqual(r['ci_method'],'normal Wald from rounded published coefficient and SE')
        self.assertFalse(r['new_raw_data_model'])
    def test_invalid_uncertainty_rejected(self):
        with self.assertRaises(ValueError): odds_ratio(1,-.1)
        with self.assertRaises(ValueError): odds_ratio(float('nan'),.1)

if __name__=='__main__': unittest.main()
