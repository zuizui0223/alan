import unittest
import pandas as pd
from scripts.analysis import alan_functional_dependence_20261115 as m

def rows(arms,**overrides):
    base={key:'x' for key in m.KEYS};base.update(study_ID='study1',plant_species='Plant one',treatment_effectiveness_metric='fruit set',effectiveness_units='proportion')
    base.update(overrides)
    return pd.DataFrame([{**base,'treatment_condition':arm,'effectiveness_value':value,'sample_size':10} for arm,value in arms.items()])
A={'day pollination':.9,'night pollination':.9,'open pollination':1.}
class DependenceTests(unittest.TestCase):
    def test_marginal_not_night_alone(self):
        p,q=m.match_panels(rows(A));r=p.iloc[0]
        self.assertAlmostEqual(r.night_loss,.1);self.assertAlmostEqual(r.night_alone,.9)
    def test_negative_deficit_retained(self):
        p,q=m.match_panels(rows({**A,'day pollination':1.2}));self.assertAlmostEqual(p.iloc[0].night_loss,-.2)
    def test_zero_open_not_divided(self):
        p,q=m.match_panels(rows({**A,'open pollination':0}));self.assertEqual(len(p),0);self.assertEqual(q['zero_open_panels'],1)
    def test_missing_arm_not_zero_filled(self):
        p,q=m.match_panels(rows({'open pollination':1,'night pollination':.9}));self.assertEqual(len(p),0)
    def test_night_missing_does_not_remove_primary_pair(self):
        p,q=m.match_panels(rows({'open pollination':1,'day pollination':.6}));self.assertEqual(len(p),1);self.assertTrue(pd.isna(p.iloc[0].night_alone))
    def test_duplicates_quarantined_not_first_row(self):
        x=rows(A);x=pd.concat([x,x.iloc[[0]]]);p,q=m.match_panels(x);self.assertEqual(len(p),0);self.assertEqual(q['ambiguous_panels'],1)
    def test_additional_treatments_never_cross_matched(self):
        x=pd.concat([rows({'open pollination':1},additional_treatment='control'),rows({'day pollination':.6},additional_treatment='light')]);p,q=m.match_panels(x);self.assertEqual(len(p),0)
    def test_endpoint_unit_never_cross_matched(self):
        x=pd.concat([rows({'open pollination':1},effectiveness_units='percent'),rows({'day pollination':.6},effectiveness_units='proportion')]);p,q=m.match_panels(x);self.assertEqual(len(p),0)
    def test_closed_reference_adjustment(self):
        p,q=m.match_panels(rows({**A,'day pollination':.8,'complete exclusion':.5}));self.assertAlmostEqual(p.iloc[0].night_loss_above_closed,.4)
    def test_invalid_closed_reference_stays_missing(self):
        p,q=m.match_panels(rows({**A,'complete exclusion':1.2}));self.assertTrue(pd.isna(p.iloc[0].night_loss_above_closed))
    def test_zero_day_night_no_false_dominance(self):
        p,q=m.match_panels(rows({**A,'day pollination':0,'night pollination':0}));self.assertTrue(pd.isna(p.iloc[0].night_share))
    def test_species_alias_and_missing_traits(self):
        p,q=m.match_panels(rows(A));tr=pd.DataFrame([{'old_name':'Plant one','plant_species':'New name','bloom_period_simple':'night'}]);r=m.join_traits(p,tr);self.assertEqual(r.iloc[0].bloom_period_simple,'night');self.assertEqual(r.iloc[0].trait_name,'New name')
    def test_conflicting_alias_raises(self):
        p,q=m.match_panels(rows(A));tr=pd.DataFrame([{'old_name':'Plant one','plant_species':'A','bloom_period_simple':'night'},{'old_name':'Plant one','plant_species':'B','bloom_period_simple':'day'}]);
        with self.assertRaises(ValueError):m.join_traits(p,tr)
    def test_study_equal_not_row_weighted(self):
        x=pd.DataFrame({'study_ID':['a']*9+['b'],'night_loss':[0]*9+[1.]});s=m.cluster_summary(x,iterations=1000,seed=7);self.assertEqual(s['studies'],2);self.assertEqual(s['median_study_median'],.5)
    def test_seed_per_fruit_not_total_reproduction(self):
        self.assertEqual(m.endpoint_subtype('seed set','seeds/fruit'),'seed_per_fruit')
        self.assertEqual(m.endpoint_subtype('seed set','seeds/plant'),'seed_per_plant')
    def test_unknown_denominator_not_claimed_per_flower(self):
        self.assertEqual(m.endpoint_subtype('fruit set','<NA>'),'fruit_unspecified')
    def test_explicit_fruit_fraction_separate_from_counts(self):
        self.assertEqual(m.endpoint_subtype('fruit set','proportion'),'fruit_proportion')
        self.assertEqual(m.endpoint_subtype('fruit set','fruit per plant'),'fruit_count')
    def test_viability_fraction_not_seed_count(self):
        self.assertEqual(m.endpoint_subtype('seed set','proportion converted (orig: percent viable seeds per fruit)'),'seed_viability_fraction')
        self.assertEqual(m.endpoint_subtype('seed set','proportion of flowers that produced seed'),'seed_formation_probability')
if __name__=='__main__':unittest.main()
