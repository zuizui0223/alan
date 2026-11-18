"""Original full Lu network source: no silent missing-light reconstruction."""
import unittest

from scripts.analysis.alan_lu2026_original_fullnetwork_sameflower_20261118 import (
  parse_source,compare_views,joint_day_night_same_plant_routes,theoretical_same_flower_fly_effectiveness_threshold,SOURCE_ROWS,COMPLETE_RECORDS
)

class TestOriginalFullNetworkSource(unittest.TestCase):
    @staticmethod
    def fixture(two_invalid=False):
        names=["time","sampling_time","treatment","block","repeat","plant",
               "all","Lepidoptera","Diptera","Hymenoptera","Coleoptera","Others"]
        out=[names]
        for j in range(COMPLETE_RECORDS):
            out.append(["1","Nocturnal","Lit" if j%2 else "Unlit",
                        str(j%5+1),"1","Plant_A",1,1,0,0,0,0])
        out.append(["1","","","1","1","Plant_A",1,1,0,0,0,0])
        if two_invalid:out[-2][1]=""
        self_count=len(out)-1
        if self_count!=SOURCE_ROWS:raise AssertionError(self_count)
        return out

    def test_single_original_unlabelled_source_record_never_forward_filled(self):
        valid,support=parse_source(self.fixture())
        self.assertEqual(len(valid),COMPLETE_RECORDS)
        self.assertEqual(support["original_rows_with_unreported_light_or_diel"],1)
        self.assertEqual(support["original_distinct_plant_labels"],1)
        self.assertTrue(support["original_unlabelled_plant_not_assumed_0_visits"])

    def test_two_missing_source_light_labels_stop_without_guessing(self):
        with self.assertRaises(ValueError):
            parse_source(self.fixture(two_invalid=True))

    def test_three_views_dont_automatically_become_one_source_cohort(self):
        valid,_=parse_source(self.fixture())
        visit=[["sampling_time","abundance","D_abundance","L_abundance",
                "H_abundance","C_abundance","O_abundance"],
               ["Nocturnal",100,1,99,0,0,0]]
        focal=[["sampling_time","all","Diptera","Lepidoptera",
                "Hymenoptera","Coleoptera","Others"],
               ["Nocturnal",10,1,9,0,0,0]]
        result=compare_views(valid,visit,focal)
        self.assertFalse(result["source_three_views_exact_same_visitation_totals"])
        self.assertTrue(result["do_not_sum_across_three_views"])

    def test_same_flower_day_night_alternative_requires_both_arms_both_periods(self):
        def row(phase,arm,plant,block,moth,fly):
            return {"phase":phase,"arm":arm,"plant":plant,"block":block,
                    "round":"t1","repeat":"1",
                    "all":moth+fly,"Lepidoptera":moth,"Diptera":fly,
                    "Hymenoptera":0,"Coleoptera":0,"Others":0}
        d=[
            row("night","dark","P","1",5,1),
            row("night","lit","P","1",1,3),
            row("day","dark","P","1",0,1),
            row("day","lit","P","1",0,4),
            row("night","dark","Q","2",4,2),
            row("night","lit","Q","2",1,2)
        ]
        r=joint_day_night_same_plant_routes(d)
        self.assertEqual(r["same_plant_and_spatial_block_has_both_arms_for_day_and_night"],1)
        self.assertEqual(r["patterns"]["night_moth_down_and_both_fly_up"]["source_plant_block_pairs"],1)
        self.assertAlmostEqual(r["positive_night_fly_visit_gain_on_those_same_plant_pairs"],2)
        self.assertAlmostEqual(r["positive_day_fly_visit_gain_on_those_same_plant_pairs"],3)
        self.assertTrue(r["same_flower_temporal_visitor_change_not_a_pollen_or_seed_rescue"])

    def test_threshold_is_counterfactual_not_observed_pollen_effectiveness(self):
        def datum(block,arm,moth,fly):
            return {"phase":"night","arm":arm,"plant":"focal","block":block,
                    "round":"t1","repeat":"r1",
                    "all":moth+fly,"Lepidoptera":moth,"Diptera":fly,
                    "Hymenoptera":0,"Coleoptera":0,"Others":0}
        d=[datum("1","dark",4,1),datum("1","lit",1,3),
           datum("2","dark",4,1),datum("2","lit",1,3)]
        x=theoretical_same_flower_fly_effectiveness_threshold(d)
        p=x["source_original_plant_specific_counterfactual_thresholds"]["focal"]
        self.assertEqual(p["observed_same_lit_dark_spatial_blocks"],2)
        self.assertAlmostEqual(p["source_moth_visits_lit_minus_dark"],-6)
        self.assertAlmostEqual(p["source_fly_visits_lit_minus_dark"],4)
        self.assertAlmostEqual(p["counterfactual_min_q_fly_per_visit_relative_to_q_moth"],1.5)
        self.assertTrue(x["not_measurement_of_fly_or_moth_pollen_efficiency"])

if __name__=="__main__":
    unittest.main()
