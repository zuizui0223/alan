"""Source-support Jaccard and algebraic matched/unmatched decomposition guards."""
import math,unittest
import pandas as pd
from scripts.analysis.alan_Giavi2021_coobserved_portfolio_decomposition_20261115 import (
    snapshot,totals,density,summarize_source
)

def minimal():
    rows=[]
    for s in range(6):
      for day in range(7 if s==0 else 6):
        pair=s*7+day
        common={'site':f's{s}','pair':pair,'plant':'common species'}
        for arm in ['C','T']:
          rows.append({**common,'treatment':arm,
            'floral_cover':10.,'all3_visits':10. if arm=='C' else 8.})
        # arm-only sample should have an actually recorded source row,
        # not be treated as biological absence in the opposite treatment.
        for arm,flower,visits in [('C',2.,1.),('T',4.,1.)]:
          rows.append({'site':f's{s}','pair':pair,
            'plant':'dark species' if arm=='C' else 'lit species',
            'treatment':arm,'floral_cover':flower,
            'all3_visits':visits})
    return pd.DataFrame(rows)

class SourceSupportTests(unittest.TestCase):
    def test_exact_support_union_and_rate_identity(self):
        v,key=snapshot(minimal())
        self.assertEqual(len(v),148)
        self.assertEqual(len(key),111)
        self.assertEqual(key.n_source_arms.eq(2).sum(),37)
        a=totals(v)
        self.assertEqual(a['both_arms_observed']['C']['n_source_plant_observations'],37)
        self.assertEqual(a['arm_only_observed']['T']['n_source_plant_observations'],37)
        self.assertAlmostEqual(density(a['both_arms_observed'])['lit_over_control_visitation_per_cover_proxy'],.8)
        self.assertAlmostEqual(density(a['arm_only_observed'])['lit_over_control_visitation_per_cover_proxy'],.5)
    def test_summary_matched_support_is_not_caused_by_missed_meadow(self):
        z=summarize_source(minimal())
        self.assertEqual(z['source_support']['both_arms_source_plant_survey_keys'],37)
        self.assertEqual(z['source_support']['exact_site_survey_plant_union'],111)
        self.assertAlmostEqual(z['source_support']['co_recorded_keys_over_union'],1/3)
        self.assertEqual(z['source_reproduction']['actual_independent_light_dark_site_pairs'],6)
        self.assertEqual(z['source_reproduction']['actual_paired_source_sampling_sessions'],37)
        self.assertTrue(z['biological_inference_gates']['arm_only_observations_HAVE_NO_matched_treatment_counterfactual'])
    def test_repeat_duplicate_fails_closed(self):
        a=minimal();a=pd.concat([a,a.iloc[[0]]],ignore_index=True)
        with self.assertRaises(ValueError):
            snapshot(a)
    def test_no_plant_record_cannot_be_relabelled_zero(self):
        g=minimal().copy()
        g.loc[0,'treatment']='T'
        with self.assertRaises(ValueError):
            snapshot(g)

if __name__=='__main__':unittest.main()
