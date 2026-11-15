"""No invented common endpoint, exact natural-log arithmetic across true studies."""
import math,unittest
from scripts.analysis.alan_P3_LiBoom_original_reproductive_pathway_accounting_20261115 import log_factors
from scripts.analysis.alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114 import (
    EXPECTED,decomposition
)
class TwoStudyAccounting(unittest.TestCase):
    def table(self):
        out={}
        for name,ref in EXPECTED.items():
            p,i,f,printed,flowers,fruits,ovules=ref
            out[name]={"fruit_set":p,
             "infestation_given_fruit":i,
             "ovule_fertilization_given_uninfested_fruit":f,
             "estimated_seeds_per_flower_printed":printed,
             "fruit_set_source_nested_flowers":flowers,
             "infestation_source_nested_fruits":fruits,
             "ovule_sample_nested_in_fruit_plant_transect_site":ovules}
        return out
    def test_two_stage_source_ratio_identity(self):
        z=log_factors({"fruit_formation":3.5572593167701863,
           "seeds_per_fruit":1.1917545541706616})
        self.assertAlmostEqual(z["product_stage_ratio"],4.239379991126886)
        self.assertAlmostEqual(sum(z["log_stage_ratio"].values()),math.log(z["product_stage_ratio"]))
    def test_primary_boom_four_arm_published_source_identity(self):
        q=decomposition(self.table())
        self.assertEqual(set(q),set(EXPECTED))
        self.assertAlmostEqual(q["Dark"]["ratio_to_dark_recomputed"],1.)
        self.assertTrue(.50<q["Green"]["ratio_to_dark_recomputed"]<.53)
        self.assertLess(q["White"]["ratio_to_dark_recomputed"],1.)
        self.assertLess(q["Red"]["ratio_to_dark_recomputed"],1.)
    def test_boom_fruit_slight_increase_can_coexist_with_seed_loss(self):
        b=decomposition(self.table())["Green"]
        f=b["algebraic_source_product_ratio_factors"]
        self.assertGreater(f["fertilized_flower_factor"],1.)
        self.assertLess(f["escape_seed_predation_factor"],1.)
        self.assertLess(f["fertilized_ovule_factor"],1.)
        self.assertAlmostEqual(math.prod(f.values()),b["ratio_to_dark_recomputed"])
    def test_log_factor_rejects_invalid_noncomparison(self):
        with self.assertRaises(ValueError):log_factors({"fruit":0,"seed":2})
        with self.assertRaises(ValueError):log_factors({"fruit":-1,"seed":2})
if __name__=="__main__":unittest.main()
