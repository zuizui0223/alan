import unittest
import pandas as pd
from scripts.analysis import alan_B_global340_paired_stages_20261115 as m
class TestGlobalDielStages(unittest.TestCase):
    def mk(self,study,stage,day,opened,year="2021",unit="number"):
        d={k:"<NA>" for k in m.EXACT}
        d.update(study_ID=study,plant_species="Test plant",year_start=year,
                 day_mean=day,open_mean=opened,stage=stage,effectiveness_units=unit)
        return d
    def test_exact_pair(self):
        p,_=m.exact_pairs(pd.DataFrame([
          self.mk("A","fruit_proportion",.3,.5),
          self.mk("A","seed_per_fruit",2,8)]))
        self.assertEqual(len(p),1)
        self.assertAlmostEqual(p[0]["delta_seed_minus_fruit"],.35)
    def test_different_years_do_not_pair(self):
        p,_=m.exact_pairs(pd.DataFrame([
          self.mk("A","fruit_proportion",.3,.5,"2020"),
          self.mk("A","seed_per_fruit",2,8,"2021")]))
        self.assertEqual(len(p),0)
    def test_ambiguous_stages_do_not_pair(self):
        p,q=m.exact_pairs(pd.DataFrame([
          self.mk("A","fruit_proportion",.3,.5),
          self.mk("A","fruit_proportion",.2,.5),
          self.mk("A","seed_per_fruit",2,8)]))
        self.assertEqual(len(p),0)
        self.assertEqual(len(q["ambiguities"]),1)
    def test_fruit_cannot_exceed_one(self):
        p,q=m.exact_pairs(pd.DataFrame([
          self.mk("A","fruit_proportion",3,5),
          self.mk("A","seed_per_fruit",2,8)]))
        self.assertEqual(len(p),0)
        self.assertEqual(q["invalid"][0]["reason"],"FRUIT_PROPORTION_GT_ONE")
    def test_light_pooled_not_strict(self):
        p,_=m.exact_pairs(pd.DataFrame([
          self.mk("Macgregor_2019","fruit_proportion",.3,.5),
          self.mk("Macgregor_2019","seed_per_fruit",2,8)]))
        self.assertEqual(len(p),1)
        self.assertTrue(p[0]["known_lighting_pooled"])
        self.assertEqual(m.summarize(p)["seed_per_fruit"]["pairs"],0)
    def test_paper_not_panel_is_sign_test_unit(self):
        p,_=m.exact_pairs(pd.DataFrame([
          self.mk("A","fruit_proportion",.3,.5,"2020"),
          self.mk("A","seed_per_fruit",2,8,"2020"),
          self.mk("A","fruit_proportion",.3,.5,"2021"),
          self.mk("A","seed_per_fruit",2,8,"2021"),
          self.mk("B","fruit_proportion",.3,.5,"2021"),
          self.mk("B","seed_per_fruit",2,8,"2021")]))
        r=m.summarize(p)["seed_per_fruit"]
        self.assertEqual((r["pairs"],r["papers"]),(3,2))
        self.assertEqual(r["posthoc_two_sided_sign_probability"],.5)
    def test_original_fixed_hash(self):
        self.assertEqual(m.SOURCE_SHA,"1d9c4c778206584130308b99ba21d902b0ccf1047118dfcf52abaec9d652e2aa")
if __name__=="__main__":unittest.main()
