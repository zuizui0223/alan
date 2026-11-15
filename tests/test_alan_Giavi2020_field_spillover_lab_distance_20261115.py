"""True published Table1/2 values synthetic regression fixtures, no plant-level fit."""
import math, unittest
from scripts.analysis.alan_Giavi2020_field_spillover_lab_distance_20261115 import (
    analyze, transform, EXPECTED
)

def fixture():
    def row(label,tag,p="0.1"):
        x,y=EXPECTED[tag]
        return [label,str(x),str(y),"0",p]
    a=[["(a)"],row("Dark (intercept)","field_epil_lit_seed_presence"),
       row("Illuminated","field_epil_lit_seed_presence"),
       ["(b)"],row("Adjacent (intercept)","field_seedpred_lit_vs_adj"),
       row("Dark","field_seedpred_dark_vs_adj"),
       row("Illuminated","field_seedpred_lit_vs_adj"),
       ["(c)"],row("Adjacent (intercept)","field_unattacked_repro_lit_vs_adj"),
       row("Dark","field_unattacked_repro_dark_vs_adj"),
       row("Illuminated","field_unattacked_repro_lit_vs_adj")]
    b=[["(a)"],row("Dark: distance from lamp","lab_visit_dark_dist_lamp"),
       row("Illuminated: distance from lamp","lab_visit_lit_dist_lamp"),
       row("Illuminated","lab_visit_lit_intercept"),
       ["(b)"],row("Dark: distance from lamp","lab_eggs_dark_dist_lamp"),
       row("Illuminated: Distance from lamp","lab_eggs_lit_dist_lamp"),
       row("Illuminated","lab_eggs_lit_intercept"),
       ["(c)"],row("Illuminated","lab_total_eggs_lit")]
    return {"article_doi":"10.1038/s41598-020-68667-y","article_xml_sha256":"source-test",
        "original_published_model_tables":[
          {"label":"Table 1","rows":a},
          {"label":"Table 2","rows":b}]}

class GiaviLabSourceTests(unittest.TestCase):
    def test_field_adjacent_actual_infestation_model_OR(self):
        z=analyze(fixture())
        x=z["original_field_stage_coefficients"]
        self.assertAlmostEqual(x["adjacent_dark_vs_far_dark_infected_fruit_odds"]["exponentiated_coefficient"],math.exp(1.90))
        self.assertAlmostEqual(x["adjacent_dark_vs_direct_lit_infected_fruit_odds"]["exponentiated_coefficient"],math.exp(1.14))
        self.assertTrue(x["unattacked_fruit_reproduction_dark_vs_adjacent_source_link"]["source_bivariate_binomial_matrix_denominator_ambiguous"])
        self.assertFalse(z["source_inference_limits"]["no_new_individuallevel_repro_model_fitted"] is False)
    def test_poisson_distance_slope_lit_positive_but_dark_not(self):
        z=analyze(fixture())["original_lab_Table2_visitation_response"]
        self.assertAlmostEqual(z["lit_distance_per_1m"]["exponentiated_coefficient"],math.exp(.37))
        self.assertAlmostEqual(z["lit_distance_per_1m"]["source_two_meter_distance_multiplier"],math.exp(.74))
        self.assertAlmostEqual(z["dark_distance_per_1m"]["exponentiated_coefficient"],math.exp(-.08))
        self.assertTrue(z["lit_distance_per_1m"]["Wald95_from_rounded_published_SE"][0]>1)
        self.assertTrue(z["lit_distance_per_1m"]["not_original_plantlevel_fitted_regression"])
    def test_eggs_and_lamp_effect_are_not_net_seed_output(self):
        z=analyze(fixture())
        q=z["original_lab_Table2_oviposition_response"]
        self.assertAlmostEqual(q["lit_distance_per_1m"]["source_two_meter_distance_multiplier"],math.exp(.48))
        self.assertAlmostEqual(q["lit_conditional_at_model_zero_distances"]["exponentiated_coefficient"],math.exp(-1.17))
        self.assertTrue(z["new_lab_original_table2_behavior"]["mechanistic_field_laboratory_consistency_not_causal_mediation"])
        self.assertEqual(z["published_lab_design"]["independent_moth_individuals"],22)
        self.assertEqual(z["published_lab_design"]["true_release_sessions"],27)
        self.assertTrue(z["published_field_design"]["adjacent_dark_and_direct_light_share_same_four_meadows"])
    def test_never_accept_changed_primary_coefficient(self):
        z=fixture()
        for v in z["original_published_model_tables"]:
            if v["label"]=="Table 2":
                for q in v["rows"]:
                    if len(q)==5 and q[0]=="Illuminated: distance from lamp":
                        q[1]="0.88"
        with self.assertRaises(ValueError):analyze(z)
    def test_wald_exp_is_correct_link_not_risk_ratio(self):
        a=transform(.37,.104,"log")
        self.assertEqual(a["model_link"],"log")
        self.assertAlmostEqual(a["exponentiated_coefficient"],math.exp(.37))
        self.assertTrue(a["Wald95_from_rounded_published_SE"][0]<math.exp(.37)<a["Wald95_from_rounded_published_SE"][1])
        with self.assertRaises(AssertionError):transform(1,.1,"identity")
if __name__=="__main__":unittest.main()
