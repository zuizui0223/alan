"""Synthetic tests for matched original 3-site Li array-vs-point source audit."""
import unittest
from scripts.analysis.alan_Li2026_original_matched_strata_bud_audit_20261114 import (
    pair_result, strata_table, by_round, source_coverage, eligible_group, missingness_tipping_bound
)

def record(sp, field, rnd, arm, yes, buds=8):
    return {"species":sp,"replicate":str(field),"round":str(rnd),"Treatment":arm,
            "pollinated":str(int(yes)),"Seeds":str(12 if yes else 0),
            "Fruit":str(2 if yes else 0),"Starting":str(buds)}

class TestMatchedOriginalEcologicalGrain(unittest.TestCase):
    def test_calculate_lit_pair_with_real_seed_vs_buds_denominators(self):
        x=[record("ML",1,1,"L",1,buds=10), record("ML",1,1,"L",0,buds=5),
           record("ML",1,1,"P",0,buds=8), record("ML",1,1,"C",0,buds=9)]
        z=pair_result(x,"L","P")
        self.assertTrue(z["comparable"])
        self.assertAlmostEqual(z["seed_presence_difference"],.5)
        self.assertAlmostEqual(z["seeds_per_starting_bud_difference"],12/15)
        self.assertAlmostEqual(z["starting_bud_mean_difference"],-.5)
        self.assertFalse(pair_result(x,"C","L")["seed_presence_difference"]>0)

    def test_missing_plant_species_arm_is_not_imputed_zero(self):
        x=[record("ML",1,1,"L",1),record("ML",1,1,"C",0),
           record("SN",1,1,"L",1)]
        z=strata_table(x,("species","replicate"),"L","C")
        self.assertEqual(z["n_comparable_strata"],1)
        self.assertEqual(z["n_unmatched"],1)
        self.assertEqual(z["n_positive"],1)
        self.assertTrue(z["not_independent_experimental_replicates"])
        self.assertFalse([a for a in z["strata"] if a["stratum"]["species"]=="SN"][0]["comparable"])

    def test_round_is_separate_composition_not_contemporaneous_census(self):
        x=[record("ML",1,1,"L",1),record("ML",1,1,"C",0),
           record("ML",1,2,"L",0),record("ML",1,2,"C",1)]
        z=by_round(x,"L","C")
        self.assertAlmostEqual(z[0]["same_round_original_plant_level"]["seed_presence_difference"],1.0)
        self.assertAlmostEqual(z[1]["same_round_original_plant_level"]["seed_presence_difference"],-1.0)
        self.assertFalse(z[2]["same_round_original_plant_level"]["comparable"])

    def test_adversarial_missingness_is_scenario_and_flips_only_at_enough_n(self):
        # Observed aggregate synthetic fixture matching author 79-plant
        # binary counts, not an asserted raw source reconstruction.
        x=[]
        for arm,n,yes in (("L",27,20),("P",25,12),("C",26,8)):
            for i in range(n):
                x.append(record("ML",1,1,arm,int(i<yes)))
        b=missingness_tipping_bound(x,"L","C",max_missing=30)
        self.assertAlmostEqual(
            b["hypothetical_11_missing_worst_case"]["minimum_possible_unadjusted_risk_difference"],
            0.18670576735092864)
        self.assertEqual(b["first_possible_sign_reversal_if_adversarial_missing"]["hypothetical_total_missing"],23)
        self.assertTrue(b["not_actual_assignment_or_missingness_model"])
        c=missingness_tipping_bound(x,"P","C",max_missing=30)
        self.assertLess(c["hypothetical_11_missing_worst_case"]["minimum_possible_unadjusted_risk_difference"],0)
        self.assertEqual(c["first_possible_sign_reversal_if_adversarial_missing"]["hypothetical_total_missing"],9)

    def test_archived_rows_differ_from_potential_design_slots(self):
        x=[record("ML",1,1,"L",1),record("SN",2,1,"L",1),record("ZC",3,1,"C",0)]
        z=source_coverage(x)
        self.assertEqual(z["total_potential_slots_described"],90)
        self.assertEqual(z["total_original_csv_rows"],3)
        self.assertTrue(z["do_not_infer_unrecorded_rows_as_failed_pollination"])

if __name__=="__main__":
    unittest.main()
