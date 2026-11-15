import csv
import importlib.util
import io
import math
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_Lu2026_iwe_access_scope_sync_20261115.py"
spec=importlib.util.spec_from_file_location("lu_crossrepo", SOURCE)
M=importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

def rows():
    csv_text="""doi,plant,endpoint,pollination_access,reported_effect,percent_change,lit_dark_ratio,endpoint_definition,source_pages,sem,sampling_variance,p_value,raw_data_obtained
"""
    data=[]
    for p in M.PLANTS:
        for e in M.ENDPOINTS:
            for a in M.ACCESSES:
                k=(p,e,a)
                if k in M.NUMERIC:
                    pct,rr=M.NUMERIC[k];label="decrease_reported"
                else:
                    pct,rr="","";label="no_significant_effect_reported"
                data.append(dict(doi=M.DOI,plant=p,endpoint=e,pollination_access=a,
                    reported_effect=label,percent_change=str(pct),lit_dark_ratio=str(rr),
                    endpoint_definition=M.DEFINITIONS[e],source_pages="14",
                    sem="",sampling_variance="",p_value="",raw_data_obtained="no"))
    return data

class TestOriginalLuAccessCorrectness(unittest.TestCase):
    def test_full_12_source_grid(self):
        v=M.validate(rows())
        self.assertEqual(v["source_combinations_count"],12)
        self.assertEqual(v["source_nonmissing_numeric_ratios"],3)
        self.assertEqual(v["nonsignificant_but_unquantified_ratios"],9)
    def test_56_percent_only_night(self):
        v=M.validate(rows())
        self.assertTrue(v["Pedicularis_seed_NIGHT_ONLY_minus_56_percent"])
        self.assertEqual(v["Pedicularis_seed_night_only_lit_dark_ratio"],.44)
        self.assertIsNone(v["Pedicularis_seed_open_lit_dark_ratio"])
        self.assertEqual(v["numeric_open_vs_night_seed_pairs"],0)
    def test_open_night_pollen_point(self):
        v=M.validate(rows())
        p=v["Pedicularis_pollen_open_vs_night"]
        self.assertAlmostEqual(p["ratio_of_response_ratios"],.7691/.462)
        self.assertAlmostEqual(p["log_ratio_of_response_ratios"],math.log(.7691/.462))
        self.assertIsNone(p["CI95"])
        self.assertFalse(p["causal_daytime_compensation_inferred"])
    def test_reject_wrong_open_seed(self):
        r=rows();one=next(x for x in r if x["plant"]=="Pedicularis kansuensis" and x["endpoint"]=="seed_set" and x["pollination_access"]=="open")
        one.update(reported_effect="decrease_reported",percent_change="-56.0",lit_dark_ratio=".44")
        with self.assertRaises(ValueError):
            M.validate(r)
    def test_reject_zero_imputation_for_nonsignificant(self):
        r=rows();one=next(x for x in r if x["plant"]=="Euphrasia regelii" and x["endpoint"]=="seed_set" and x["pollination_access"]=="open")
        one["lit_dark_ratio"]="1"
        with self.assertRaises(ValueError):
            M.validate(r)
    def test_no_meta_slope(self):
        v=M.validate(rows())
        self.assertIsNone(v["global_ALAN_night_dependence_beta_V"])
        self.assertEqual(v["independent_light_by_access_reproduction_meta_effects"],0)
    def test_reject_wrong_denominator(self):
        r=rows();r[0]["endpoint_definition"]="seeds_per_plant"
        with self.assertRaises(ValueError):
            M.validate(r)

if __name__=="__main__":
    unittest.main()
