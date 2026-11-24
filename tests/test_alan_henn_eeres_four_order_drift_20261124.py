"""Source vocabulary and physical-block estimand invariants for two programmes."""
import importlib.util
import math
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analysis/alan_20261124_henn_eeres_taxon_matched_drift.py"
spec=importlib.util.spec_from_file_location("matched4",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestMatchedDriftComposition(unittest.TestCase):
    def test_all_source_groups_unique_and_accounted_for(self):
        groups=mod.TAXON_GROUPS
        self.assertEqual({k:len(v) for k,v in groups.items()},
                         {"E_T":11,"D_C":18,"Plecoptera_excluded":3,"NonInsect_excluded":6})
        names=[name for s in groups.values() for name in s]
        self.assertEqual(len(names),38)
        self.assertEqual(len(set(names)),38)
        mapped=mod.source_taxon_mapping(names)
        self.assertEqual(len(mapped),38)
        self.assertEqual(mapped["Baetis"],"E_T")
        self.assertEqual(mapped["Chironomini"],"D_C")
        self.assertEqual(mapped["Gammarus"],"NonInsect_excluded")
        self.assertEqual(mapped["Perlodes"],"Plecoptera_excluded")
        with self.assertRaises(ValueError):
            mod.source_taxon_mapping(names[:-1]+["unknown_taxon"])
        with self.assertRaises(ValueError):
            mod.source_taxon_mapping(names[:-1]+[names[0]])

    def test_five_independent_henn_rivers_source_derived_order_replay(self):
        a=mod.henn_study()
        self.assertEqual(a["n_independent_physical_units"],5)
        self.assertEqual(a["table_total_discrepancy_lit_count"],2)
        self.assertAlmostEqual(a["four_order_composition_RRR"],1.133273435160227,places=12)
        self.assertAlmostEqual(a["delta_SE_shared_physical_units"],.07847123889548929,places=12)
        self.assertAlmostEqual(a["descriptive_t95"][0],.91141149,places=7)
        self.assertAlmostEqual(a["descriptive_t95"][1],1.40914252,places=7)
        self.assertAlmostEqual(a["alternate_equal_stream_geometric_RRR"],1.4024150760152,places=12)
        self.assertEqual(a["streams_with_four_order_RRR_above_one"],3)
        self.assertEqual(len(a["block_marginals"]),5)
        self.assertEqual(a["block_marginals"][-1],{"A_lit":560,"B_lit":174,"A_dark":837,"B_dark":252})

    def test_same_study_not_pseudoreplicated(self):
        a=mod.ratio_from_block_arms([(10.,5.,5.,5.) for _ in range(6)],"mock")
        self.assertEqual(a["n_independent_physical_units"],6)
        self.assertAlmostEqual(a["four_order_composition_RRR"],2.0)
        self.assertEqual(a["delta_SE_shared_physical_units"],0)
        with self.assertRaises(ValueError):
            mod.ratio_from_block_arms([(10.,5.,5.,5.) for _ in range(8)],"fake24")
        with self.assertRaises(ValueError):
            mod.ratio_from_block_arms([(0.,0.,0.,0.) for _ in range(6)],"invalid")


if __name__=="__main__":
    unittest.main()
