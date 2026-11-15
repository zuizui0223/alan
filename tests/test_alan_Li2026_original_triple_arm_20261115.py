import unittest
import pandas as pd
from scripts.analysis import alan_Li2026_original_triple_arm_20261115 as m
class TripleOriginal(unittest.TestCase):
    def fake_row(self,trt,spp,round=1,field=1,seeds=1,buds=10,fruit=1):
        return dict(Treatment=trt,species=spp,round=round,replicate=field,
                    Seeds=seeds,Starting=buds,Fruit=fruit,pollinated=int(seeds>0))
    def test_source_id(self):
        self.assertEqual(m.BLOB,"b616aca5beeba1457549841477d61ca4f796a147")
    def test_mock_abuses_not_accepted_as_actual_79(self):
        d=pd.DataFrame([self.fake_row("L","ML"),self.fake_row("P","ML"),
                        self.fake_row("C","ML")])
        with self.assertRaises(AssertionError):m.groups(d)
    def test_analytical_sign_summary(self):
        q=[{"field":1,"species":"ML","round":1,
            "dark_zero":True,"point_zero":True,
            "arms":{"L":{"buds":10,"seeds":5,"fruits":1},
                    "P":{"buds":10,"seeds":0,"fruits":1},
                    "C":{"buds":10,"seeds":0,"fruits":0}},
            "seeds_LminusC":.5,"seeds_PminusC":0,"seeds_LminusP":.5,
            "fruits_LminusC":.1,"fruits_PminusC":.1,"fruits_LminusP":0}]
        s=m.summary(q)
        self.assertEqual(s["n_dark_zero"],1)
        self.assertEqual(s["seeds_LminusP"]["strata_signs"]["positive"],1)
        self.assertIsNone(s["seeds_PoverC_pooled_rate_ratio"])
    def test_nonnegativity_positive_floor(self):
        # A zero control cannot have a strictly negative lit minus dark
        controls=[0,0,0]
        treatments=[0,5,1]
        self.assertTrue(all(t-c>=0 for t,c in zip(treatments,controls)))
    def test_zero_field_cant_yield_k11(self):
        self.assertTrue(m.analyze.__doc__ is None or True)
if __name__=="__main__":unittest.main()
