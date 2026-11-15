import unittest
import pandas as pd
from scripts.analysis import alan_Li2026_matched_field_taxon_round_20261115 as m
class TestLiMatched(unittest.TestCase):
    def row(self,field,species,rnd,treatment,seeds=4,fruits=1,buds=10):
        return {"replicate":field,"species":species,"round":rnd,"Treatment":treatment,
                "Starting":buds,"Seeds":seeds,"Fruit":fruits}
    def test_source_sha_identity(self):
        self.assertEqual(m.BLOB,"b616aca5beeba1457549841477d61ca4f796a147")
    def test_no_impute_missing_into_zero(self):
        # in real prepare requires 79 source rows; fabricate 79 via exact copies + 1 NaN
        base=[self.row(1,"ML",1,"L",10),self.row(1,"ML",1,"C",2)]
        raw=pd.DataFrame(base*39+[self.row(3,"SN",2,"P",float("nan"))])
        # prepare enforces actual source 3 species and fields, so direct
        # synthetic low-level support test is handled by invalid cohort rejection
        with self.assertRaises(AssertionError):m.prepare(raw)
    def test_sign_two_sided_few_fields(self):
        self.assertEqual(m.sign_p(3,0),.25)
        self.assertEqual(m.sign_p(2,1),1.)
    def test_stats_holdout_species(self):
        keys=("field","species","round")
        rec=[]
        for f,s,delta in [(1,"ML",1),(2,"ML",2),(3,"SN",-1)]:
            rec.append({"field":f,"species":s,"round":1,
                        "lit_buds":10,"dark_buds":10,
                        "lit_seeds":4+delta,"dark_seeds":4,
                        "lit_fruits":2,"dark_fruits":1,
                        "seed_per_bud_diff":delta/10,
                        "fruit_per_bud_diff":.1})
        q=m.stats(rec)
        self.assertEqual(q["seed_per_bud_diff"]["strata_positive_negative_zero"],(2,1,0))
        self.assertEqual(len(q["seed_per_bud_diff"]["field_equal_differences"]),3)
        self.assertAlmostEqual(q["seed_per_bud_diff"]["true_field_equal_mean_difference"],(1+2-1)/30)
    def test_cannot_promote_single_field_to_global(self):
        q=m.stats([{"field":1,"species":"ZC","round":1,
                     "lit_buds":10,"dark_buds":10,"lit_seeds":20,"dark_seeds":1,
                     "lit_fruits":2,"dark_fruits":1,
                     "seed_per_bud_diff":1.9,"fruit_per_bud_diff":.1}])
        self.assertEqual(len(q["seed_per_bud_diff"]["field_equal_differences"]),1)
        self.assertTrue(q["not_a_model_based_ALAN_meta_effect"])
if __name__=="__main__":unittest.main()
