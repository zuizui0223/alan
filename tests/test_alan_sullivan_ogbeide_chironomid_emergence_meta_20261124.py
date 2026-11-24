"""Physical-unit taxon-matched emergence meta, source-independent offline unit tests."""
import importlib.util,copy,math,json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analysis/alan_sullivan_ogbeide_chironomid_emergence_meta_20261124.py"
spec=importlib.util.spec_from_file_location("matched_emergence",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def synthetic_sullivan():
    # Test-only fictional allocation, with authentic *total* constraints;
    # this is NOT the actual Sullivan per-block empirical data.
    counts={
        "DRK":[120,118,118],
        "LEDI":[400,400,362],
        "LEDH":[500,400,395],
        "HPSI":[500,500,457],
        "HPSH":[200,230,260],
    }
    blocks=[dict(group=str(i+1),light=arm,
                 all_aquatic_capture=count+20,
                 explicit_chironomid_recorded_capture=count,
                 four_original_dates=4)
            for i in range(3) for arm,x in counts.items() for count in [x[i]]]
    return {
        "status":"SOURCE_LABEL_RETRIEVED_NO_META_EFFECT",
        "original_2013_author_TSV_sha256":m.sullivan.DIGEST,
        "chironomid_source_labels":["Chironomidae"],
        "original_stations":15,
        "independent_habitat_blocks":3,
        "chironomid_count_by_treatment":{tr:sum(c) for tr,c in counts.items()},
        "source_block_treatment_counts":blocks,
    }

class TestChironomidEmergenceMeta(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=synthetic_sullivan()
        cls.flumes,cls.sha,cls.original_week4=m.ogbeide.source_verify_and_read()

    def test_source_three_blocks_and_four_shared_dark_lamp_comparisons(self):
        x,variants,units=m.independent_habitat_block_estimates(self.source)
        self.assertEqual(x["true_block_n"],3)
        self.assertEqual(len(x["three_block_lnRR"]),3)
        self.assertEqual(len(variants),7)
        self.assertEqual(len(units),15)
        self.assertTrue(all(q>0 for q in units.values()))
        for g in ("1","2","3"):
            logs=[math.log(units[f"{g}/{lamp}"]/units[f"{g}/DRK"])
                  for lamp in ("LEDI","LEDH","HPSI","HPSH")]
            self.assertAlmostEqual(x["three_block_lnRR"][int(g)-1],
                                   sum(logs)/4,places=12)

    def test_source_missing_taxon_and_dark_control_stop(self):
        s=copy.deepcopy(self.source)
        s["chironomid_source_labels"]=[]
        with self.assertRaises(ValueError):m.independent_habitat_block_estimates(s)
        s=copy.deepcopy(self.source)
        s["source_block_treatment_counts"][0]["explicit_chironomid_recorded_capture"]=0
        with self.assertRaises(ValueError):m.independent_habitat_block_estimates(s)
        s=copy.deepcopy(self.source)
        s["original_stations"]=30
        with self.assertRaises(ValueError):m.independent_habitat_block_estimates(s)

    def test_one_ogbeide_study_two_weeks_4_lit_4_dark(self):
        x,variants=m.independent_flume_estimates(self.flumes)
        self.assertEqual((x["original_flume_n_lit"],x["original_flume_n_dark"]),(4,4))
        self.assertEqual(len(x["flume_mean_values_by_treatment"]["alan"]),4)
        self.assertEqual(len(x["flume_mean_values_by_treatment"]["control"]),4)
        self.assertAlmostEqual(x["RR"],
                ((83.46875+22.78125)/2)/((94.48958333333333+33.19791666666667)/2),
                places=12)
        self.assertAlmostEqual(variants["Week 4"]["RR"],.6862252902416064,places=12)
        self.assertTrue(0<x["SE"]<1)

    def test_frozen_source_receipt_stops_unexpected_reanalysis(self):
        file=ROOT/"database/alan_experiments/arthropod_expansion_20261122/chironomid_emergence_20261124/chironomid_emergence_k2_frozen_20261124.json"
        frozen=json.loads(file.read_text(encoding="utf-8"))
        self.assertEqual(len(frozen["programme_estimates"]),2)
        self.assertEqual(frozen["exploratory_meta"]["k_independent_programmes"],2)
        self.assertAlmostEqual(frozen["programme_estimates"][0]["RR"],2.662188561008355,places=12)
        self.assertAlmostEqual(frozen["programme_estimates"][1]["RR"],.8321096426823299,places=12)
        self.assertAlmostEqual(frozen["exploratory_meta"]["pooled_RR"],1.4062321732688763,places=12)
        self.assertLess(frozen["exploratory_meta"]["t1_95"][0],1)
        self.assertGreater(frozen["exploratory_meta"]["t1_95"][1],1)
        synthetic=m.compute(self.source,self.flumes)
        with self.assertRaises(ValueError):
            m.verify_frozen(synthetic,frozen)
        # The exact source replay in Actions must supply the full empirical
        # result, this fixture checks only the validator's identity contract.
        copy_result=copy.deepcopy(synthetic)
        copy_result["Sullivan_color_menu_frozen_primary"].update({
            "RR":frozen["programme_estimates"][0]["RR"],
            "SE":frozen["programme_estimates"][0]["SE_logRR"],
        })
        copy_result["Ogbeide_no_crayfish_eight_flumes_frozen_primary"].update({
            "RR":frozen["programme_estimates"][1]["RR"],
            "SE":frozen["programme_estimates"][1]["SE_logRR"],
        })
        copy_result["exploratory_k2_random_effects_meta"].update({
            "REML_tau2_closed_form_k2":frozen["exploratory_meta"]["REML_tau2"],
            "pooled_RR":frozen["exploratory_meta"]["pooled_RR"],
            "pooled_logRR_SE_mKH":frozen["exploratory_meta"]["SE_logRR_mKH"],
            "modified_HK_t1_95":frozen["exploratory_meta"]["t1_95"],
        })
        ref=frozen["within_programme_physical_unit_diagnostics"]
        copy_result["leave_one_physical_unit_out_point_sensitivities"].update({
            "Sullivan_n3_group_ratios":list(ref["Sullivan_3_group_RRs"]),
            "Sullivan_leave_one_group_out_RR_point_only":{
                f"group_{i+1}":v for i,v in enumerate(ref["Sullivan_leave_one_group_out_RRs"])},
            "Ogbeide_eight_flume_leave_one_out_RR_point_only":{
                **{f"alan_index_{i+1}":v for i,v in enumerate(ref["Ogbeide_8_flume_leave_one_out_RRs"][:4])},
                **{f"control_index_{i+1}":v for i,v in enumerate(ref["Ogbeide_8_flume_leave_one_out_RRs"][4:])}},
            "Ogbeide_Week1_RR":ref["Ogbeide_Week1_RR"],
            "Ogbeide_Week4_RR":ref["Ogbeide_Week4_RR"],
            "between_programme_ratio_of_LitDark_RR_point_only":ref["between_programme_ratio_of_RRs_point"],
        })
        self.assertTrue(m.verify_frozen(copy_result,frozen))
        copy_result["leave_one_physical_unit_out_point_sensitivities"]["Sullivan_n3_group_ratios"][0]+=0.2
        with self.assertRaises(ValueError):
            m.verify_frozen(copy_result,frozen)

    def test_k2_reml_and_mkh_df1_bounds_not_overconfident(self):
        s,v,_=m.independent_habitat_block_estimates(self.source)
        o,w=m.independent_flume_estimates(self.flumes)
        z=m.k2_meta(s,o)
        self.assertEqual(z["k_independent_programmes"],2)
        self.assertEqual(z["df_for_meta_interval"],1)
        self.assertEqual(len(z["modified_HK_t1_95"]),2)
        self.assertTrue(z["modified_HK_t1_95"][0]<z["pooled_RR"]<z["modified_HK_t1_95"][1])
        self.assertTrue(z["not_confirmatory"])
        self.assertGreaterEqual(z["REML_tau2_closed_form_k2"],0)
        full=m.compute(self.source,self.flumes)
        self.assertEqual(full["new_exploratory_endpoint_meta_lane"],1)
        self.assertFalse(full["added_to_existing_108_effect_catalogue"])
        self.assertFalse(full["changed_existing_four_meta_effects"])

if __name__=="__main__":
    unittest.main()
