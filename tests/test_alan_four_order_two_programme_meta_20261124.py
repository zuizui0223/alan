import copy
import importlib.util
import json
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124"
FILE=HERE/"two_programme_four_order_meta.py"
spec=importlib.util.spec_from_file_location("four_order_meta", FILE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestTwoProgrammeMeta(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((HERE/"two_programme_four_order_source_estimates_20261124.json").read_text(encoding="utf-8"))
        cls.out=mod.synthesize(cls.data)

    def test_k2_closed_form_REML_and_mKH_df1_not_normally_narrow(self):
        x=self.out
        self.assertEqual(x["k_independent_programmes"],2)
        self.assertEqual(x["t_reference_df"],1)
        self.assertAlmostEqual(x["reml_tau2"],.001359904329136968,places=13)
        self.assertAlmostEqual(x["pooled_lnRRR"],.1955635325401213,places=12)
        self.assertAlmostEqual(x["pooled_RRR"],1.2159960467857396,places=12)
        self.assertAlmostEqual(x["pooled_mkh_95_interval"][0],.6398196602507326,places=12)
        self.assertAlmostEqual(x["pooled_mkh_95_interval"][1],2.3110361835694366,places=12)
        self.assertAlmostEqual(x["mkh_variance_factor"],1.0,places=12)
        self.assertTrue(x["pooled_mkh_95_interval"][0]<1<x["pooled_mkh_95_interval"][1])

    def test_original_source_lit_two_count_sensitivity(self):
        lo,hi=self.out["two_count_Henn_reallocation_pooled_point_range"]
        self.assertAlmostEqual(lo,1.2138293806370006,places=12)
        self.assertAlmostEqual(hi,1.2200548821901074,places=12)
        self.assertTrue(lo<self.out["pooled_RRR"]<hi)
        self.assertTrue(self.out["input_but_not_added_to_main_effect_catalogue"])

    def test_reject_duplicate_study_or_wrong_taxa_or_negative_se(self):
        z=copy.deepcopy(self.data)
        z["studies"][1]["programme_id"]=z["studies"][0]["programme_id"]
        with self.assertRaises(ValueError):mod.synthesize(z)
        z=copy.deepcopy(self.data)
        z["denominator_orders"]=["Diptera","Plecoptera"]
        with self.assertRaises(ValueError):mod.synthesize(z)
        z=copy.deepcopy(self.data)
        z["studies"][0]["se_paired_physical_unit_delta"]=-.01
        with self.assertRaises(ValueError):mod.synthesize(z)

    def test_replay_provenance_checks(self):
        z=self.data
        replay={"status":"SOURCE_AUTHENTICATED_EXPLORATORY_TWO_PROGRAMME_COMPARISON_NOT_ACCEPTED_META",
                "meta_k_added":0,"mkh_or_reml_pooled_not_fitted":True,
                "henn":{"four_order_log_composition_RRR":z["studies"][0]["ln_rrr"],
                        "delta_SE_shared_physical_units":z["studies"][0]["se_paired_physical_unit_delta"],
                        "four_order_composition_RRR":z["studies"][0]["rrr"],
                        "n_independent_physical_units":5},
                "eeres":{"four_order_log_composition_RRR":z["studies"][1]["ln_rrr"],
                        "delta_SE_shared_physical_units":z["studies"][1]["se_paired_physical_unit_delta"],
                        "four_order_composition_RRR":z["studies"][1]["rrr"],
                        "n_independent_physical_units":6,
                        "source_sha256":z["eeres_dryad_workbook_sha256"]}}
        self.assertTrue(mod.review_source_lock(z,replay))
        replay["eeres"]["delta_SE_shared_physical_units"]+=.01
        with self.assertRaises(ValueError):mod.review_source_lock(z,replay)


if __name__=="__main__":
    unittest.main()
