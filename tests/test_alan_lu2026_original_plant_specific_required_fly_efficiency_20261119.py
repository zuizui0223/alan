import copy
import math
import unittest
from scripts.analysis.alan_lu2026_original_plant_specific_required_fly_efficiency_20261119 import (
    source_plant_q_thresholds, num, SCHEMA, AUTHOR_SHA
)


def plant(name,treatment,moth,fly,block=1):
    return {
        "plant":name,"sampling_time":"Nocturnal","treatment":treatment,
        "block":block,"Lepidoptera":moth,"Diptera":fly,"all":moth+fly,
    }


def visit(treatment,moth,fly):
    return {"sampling_time":"Nocturnal","treatment":treatment,
            "block":1,"L_abundance":moth,"D_abundance":fly}


class OriginalLuPlantQThresholdTests(unittest.TestCase):
    def original_shape(self):
        plants=[
          plant("A","Unlit",10,1),plant("A","Lit",4,4),  # -6 moth,+3 fly, q>=2
          plant("B","Unlit",2,2),plant("B","Lit",1,8),   # -1 moth,+6 fly, q>=1/6
          plant("C","Unlit",2,0),plant("C","Lit",1,0),   # no fly gain; cannot compensate
          plant("D","Unlit",2,2), # lit undocumented: NOT a visit count zero
        ]
        visits=[visit("Unlit",16,5),visit("Lit",6,12)]
        return plants,visits

    def test_source_plant_threshold_no_pollen_or_seed_inference(self):
        plants,visits=self.original_shape()
        a=source_plant_q_thresholds(plants,visits,strict=False)
        self.assertEqual(a["schema"],SCHEMA)
        self.assertEqual(a["source_original_focal_plant_species"],4)
        self.assertEqual(a["original_documented_nocturnal_moth_down_fly_up_plant_n"],2)
        self.assertEqual(a["source_field_blocks"],5)
        self.assertEqual(a["community_visit_event_losses_moth"],10)
        self.assertEqual(a["community_visit_event_gains_fly"],7)
        self.assertAlmostEqual(a["community_conditional_qfly_qmoth_threshold_if_all_visits_act_on_same_equivalent_plants"],10/7)
        q=a["focal_six_plant_required_efficiency_ratio_descriptive_summary"]
        self.assertEqual(q["n_plant_taxa_with_both_source_changes"],2)
        self.assertAlmostEqual(q["minimum"],1/6)
        self.assertAlmostEqual(q["median"],(2+1/6)/2)
        self.assertEqual(q["maximum"],2)
        self.assertEqual(q["count_where_equal_q_visits_would_suffice"],1)
        self.assertEqual(q["count_where_fly_q_must_exceed_moth_q"],1)
        self.assertTrue(q["theoretical_only_NOT_compensated_plant_proportion"])
        self.assertTrue(a["guard"]["no_causal_pollination_rescue_fraction_reported"])
        self.assertTrue(a["guard"]["six_plant_thresholds_one_five_block_experiment_NOT_six_independent_light_studies"])
        self.assertNotIn('"A"',str(a))
        self.assertNotIn('"D"',str(a))

    def test_negative_or_absent_fly_replacement_not_forced_into_compensated_fraction(self):
        plants,visits=self.original_shape()
        for r in visits:
            if r["treatment"]=="Lit":
                r["D_abundance"]=0
        a=source_plant_q_thresholds(plants,visits,strict=False)
        self.assertIsNone(a["community_conditional_qfly_qmoth_threshold_if_all_visits_act_on_same_equivalent_plants"])
        self.assertEqual(a["original_documented_nocturnal_moth_down_fly_up_plant_n"],2)
        self.assertTrue(a["guard"]["cross_plant_community_q_threshold_not_biological_same_plant_pollination"])

    def test_missing_original_visits_do_not_become_zero(self):
        plants,visits=self.original_shape()
        plants=[r for r in plants if r["plant"]!="B"]
        a=source_plant_q_thresholds(plants,visits,strict=False)
        self.assertEqual(a["source_original_focal_plant_species"],3)
        self.assertEqual(a["original_documented_nocturnal_moth_down_fly_up_plant_n"],1)
        self.assertTrue(a["guard"]["not_equal_per_flower_visitation_effort_or_flower_availability"])

    def test_source_2026_author_SHA_and_strict_data_volume_guards(self):
        self.assertEqual(AUTHOR_SHA,"ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f")
        plants,visits=self.original_shape()
        with self.assertRaisesRegex(ValueError,"Original 1038"):
            source_plant_q_thresholds(plants,visits,strict=True)
        bad=copy.deepcopy(plants)
        bad[0]["Lepidoptera"]=-1
        with self.assertRaisesRegex(ValueError,"Invalid original visitor"):
            source_plant_q_thresholds(bad,visits,strict=False)
        bad=copy.deepcopy(plants)
        bad[0]["block"]=6
        with self.assertRaisesRegex(ValueError,"Five original paired"):
            source_plant_q_thresholds(bad,visits,strict=False)
        with self.assertRaises(ValueError):
            num(None)
        with self.assertRaises(ValueError):
            num(-1.2)


if __name__=="__main__":
    unittest.main()
