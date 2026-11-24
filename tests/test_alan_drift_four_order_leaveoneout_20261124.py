"""Check the physical-programme design and Henn order-level omission identities."""
from pathlib import Path
import importlib.util
import copy
import json
import math
import unittest

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"scripts/analysis/alan_drift_four_order_leaveoneout_20261124.py"
spec=importlib.util.spec_from_file_location("drift_leaveoneout",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class TestOrderSensitivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows,cls.main=m.Henn()

    def test_source_derived_original_four_order_baseline(self):
        self.assertEqual(len(self.rows),5)
        a=m.base.ratio_from_block_arms(m.build_arms(self.rows,None),"five_rivers")
        self.assertAlmostEqual(a["four_order_composition_RRR"],1.133273435160227,places=12)
        self.assertAlmostEqual(a["delta_SE_shared_physical_units"],.07847123889548929,places=12)

    def test_no_order_deliberately_relabelled(self):
        expected={
            "Ephemeroptera":.881652344338911,
            "Trichoptera":1.1578399913646744,
            "Diptera":.9480049383673648,
            "Coleoptera":1.2311910057085462,
        }
        for removed,want in expected.items():
            a=m.base.ratio_from_block_arms(m.build_arms(self.rows,removed),"five_rivers")
            self.assertEqual(a["n_independent_physical_units"],5)
            self.assertAlmostEqual(a["four_order_composition_RRR"],want,places=12)

    def test_exact_reml_plus_mkh_reproduces_prior_k2_main(self):
        a=m.base.ratio_from_block_arms(m.build_arms(self.rows,None),"five_rivers")
        z={"four_order_log_composition_RRR":.23181424616015578,
           "delta_SE_shared_physical_units":.05008183522476254}
        x=m.mKH_two(a,z)
        self.assertEqual(x["k_independent_programmes"],2)
        self.assertAlmostEqual(x["pooled_RRR"],1.2159960467857396,places=12)
        self.assertAlmostEqual(x["mKH_t1_95"][0],.6398196602507326,places=12)
        self.assertAlmostEqual(x["mKH_t1_95"][1],2.3110361835694366,places=12)

    def test_source_result_lock_rejects_any_tamper(self):
        frozen=json.loads((ROOT/"database/alan_experiments/arthropod_expansion_20261122/stream_drift_20261124/four_order_omission_frozen_20261124.json").read_text(encoding="utf-8"))
        pseudo={
            "source_eeres_sha256":frozen["eeres_source_sha256"],
            "comparisons":[{
                "deleted_order":r["deleted_order"],
                "henn":{"RRR":r["henn_RRR"]},
                "eeres":{"RRR":r["eeres_RRR"]},
                "exploratory_two_programme_mKH":{
                    "pooled_RRR":r["two_programme_pooled_RRR"],
                    "mKH_t1_95":r["two_programme_mKH_t1_95"],
                }
            } for r in frozen["rows"]],
        }
        self.assertTrue(m.verify_frozen(pseudo,frozen))
        broken=copy.deepcopy(pseudo)
        broken["comparisons"][3]["eeres"]["RRR"]+=.1
        with self.assertRaises(ValueError):
            m.verify_frozen(broken,frozen)
        broken=copy.deepcopy(pseudo)
        broken["comparisons"][2]["exploratory_two_programme_mKH"]["mKH_t1_95"][0]+=.1
        with self.assertRaises(ValueError):
            m.verify_frozen(broken,frozen)

    def test_no_unverified_taxa_or_physical_pseudoreplication(self):
        with self.assertRaises(ValueError):
            m.build_arms(self.rows,"Megaloptera")
        for order in m.ORDERS:
            arms=m.build_arms(self.rows,order)
            self.assertEqual(len(arms),5)
            self.assertTrue(all(all(v>=0 for v in r) for r in arms))


if __name__=="__main__":
    unittest.main()
