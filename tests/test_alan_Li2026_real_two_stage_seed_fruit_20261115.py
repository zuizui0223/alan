"""Original-reproductive-stage mathematics: tests do not invent source outcomes."""
import math,unittest
from scripts.analysis.alan_Li2026_real_two_stage_seed_fruit_20261115 import records,two_stage,analyze

def sample(seeds,fruit,bud,arm="L"):
    return {"pollinated":str(int(seeds>0)),"Seeds":str(seeds),
        "Fruit":str(fruit),"Starting":str(bud),
        "Treatment":arm,"species":"SN","replicate":"1","round":"1"}

class LiTrueReproductiveStageTests(unittest.TestCase):
    def test_exact_hurdle_and_fruit_components(self):
        x=records([sample(0,0,3),sample(10,2,5),sample(5,1,3)])
        self.assertEqual(x["seed_positive_plants"],2)
        self.assertAlmostEqual(x["probability_any_seeds"],2/3)
        self.assertAlmostEqual(x["mean_seeds_conditional_on_any_seeds"],7.5)
        self.assertAlmostEqual(x["mean_fruits_per_initial_flower_bud"],3/11)
        self.assertAlmostEqual(x["mean_seeds_per_fruit"],5.)
        self.assertAlmostEqual(x["mean_seeds_per_initial_flower_bud"],15/11)
        self.assertAlmostEqual(x["hurdle_identity_absolute_error"],0.)
        self.assertAlmostEqual(x["fruit_identity_absolute_error"],0.)
    def test_exact_between_arm_ratios_not_mediator_effects(self):
        a=records([sample(0,0,5),sample(8,2,5)])
        b=records([sample(2,1,5),sample(2,1,5)])
        q=two_stage(a,b)
        self.assertAlmostEqual(q["seed_positive_probability_ratio"],.5)
        self.assertAlmostEqual(q["conditional_seeds_given_positive_ratio"],4.)
        self.assertAlmostEqual(q["observed_seeds_per_known_plant_ratio"],2.)
        self.assertAlmostEqual(q["observed_seeds_per_initial_bud_ratio"],2.)
        self.assertAlmostEqual(q["probability_and_conditional_seed_identity_error"],0.)
        self.assertAlmostEqual(q["fruit_and_seeds_per_fruit_identity_error"],0.)
        self.assertTrue(q["noncausal_arithmetic_not_statistical_mediation"])
    def test_does_not_convert_unknown_seed_to_failure(self):
        q=[sample(1,1,1),{**sample(0,0,5),"pollinated":"NA","Seeds":"NA","Fruit":"NA"}]
        x=records(q)
        self.assertEqual((x["source_plant_rows"],x["known_seed_outcome_plants"],
                          x["missing_reproductive_outcomes"]),(2,1,1))
    def test_source_count_and_fruit_quality_fail_closed(self):
        with self.assertRaises(ValueError):
            records([sample(10,3,2)])
        with self.assertRaises(ValueError):
            records([sample(7,0,9)])
        with self.assertRaises(ValueError):
            analyze([sample(0,0,9)])
if __name__=="__main__":unittest.main()
