"""Never treat taxonomic composition as extra experimental independence."""
import copy
import importlib.util
import math
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"scripts/analysis/alan_sullivan2019_chiro_versus_other_emergence_20261124.py"
spec=importlib.util.spec_from_file_location("sullivan_test_comp",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def synthetic_source():
    counts={
        "DRK":[120,118,118],
        "LEDI":[400,400,362],
        "LEDH":[500,400,395],
        "HPSI":[500,500,457],
        "HPSH":[200,230,260],
    }
    # Sum non-Chiro = 412, including dark51 + lit361.
    other={
        "DRK":[17,17,17],
        "LEDI":[20,15,15],
        "LEDH":[25,20,22],
        "HPSI":[36,36,36],
        "HPSH":[45,45,46],
    }
    rows=[]
    for i in range(3):
        for arm in ("DRK",)+m.LIGHTS:
            k=int(counts[arm][i])
            non=int(other[arm][i])
            rows.append({
                "group":str(i+1),"light":arm,
                "all_aquatic_capture":k+non,
                "explicit_chironomid_recorded_capture":k,
                "four_original_dates":4,
            })
    return {
        "status":"SOURCE_LABEL_RETRIEVED_NO_META_EFFECT",
        "original_2013_author_TSV_sha256":m.mod.DIGEST,
        "chironomid_source_labels":["Chironomidae"],
        "independent_habitat_blocks":3,
        "original_stations":15,
        "chironomid_count_by_treatment":{k:sum(a) for k,a in counts.items()},
        "source_block_treatment_counts":rows
    }

class TestSourceComposition(unittest.TestCase):
    def test_source_unit_and_non_chiro_denominator_are_not_extra_studies(self):
        x=m.evaluate(synthetic_source())
        self.assertEqual(x["true_independent_habitat_block_n"],3)
        self.assertEqual(x["source_all_2013_aquatic_emergence_count"],5372)
        self.assertEqual(x["source_chironomid_count"],4960)
        self.assertEqual(x["dark_capture_source_counts"]["chiro"],356)
        self.assertEqual(x["four_lamp_aggregate_source_counts"]["chiro"],4604)
        self.assertEqual(x["dark_capture_source_counts"]["other"],51)
        self.assertEqual(x["four_lamp_aggregate_source_counts"]["other"],361)
        self.assertEqual(len(x["three_true_physical_block_diagnostics"]),3)
        self.assertTrue(x["composition_diagnostic_is_a_different_estimand_from_2013_emergence_meta"])
        self.assertEqual(x["new_meta_independent_k_added"],0)
    def test_source_posthoc_t2_and_all_signflip_symmetry(self):
        x=m.evaluate(synthetic_source())
        self.assertTrue(x["block_t2_descriptive_RRR_95"][0]<x["block_equal_weight_RRR"]<x["block_t2_descriptive_RRR_95"][1])
        self.assertAlmostEqual(
          x["pooled_count_ratio_of_RRs_descriptive"],
          (4604/356)/(361/51),places=12)
        self.assertGreaterEqual(x["three_block_signflip_reference_p"],0.25)
        self.assertLessEqual(x["three_block_signflip_reference_p"],1)
    def test_organism_origin_or_zeros_cannot_be_made_up(self):
        z=synthetic_source()
        z["source_block_treatment_counts"][0]["explicit_chironomid_recorded_capture"]=0
        with self.assertRaises(ValueError):m.evaluate(z)
        z=synthetic_source()
        z["source_block_treatment_counts"][0]["all_aquatic_capture"]=0
        with self.assertRaises(ValueError):m.evaluate(z)
        z=synthetic_source()
        z["source_block_treatment_counts"][0]["four_original_dates"]=3
        with self.assertRaises(ValueError):m.evaluate(z)
        z=synthetic_source()
        z["chironomid_source_labels"]=["UnidentifiedDiptera"]
        with self.assertRaises(ValueError):m.evaluate(z)

if __name__=="__main__":
    unittest.main()
