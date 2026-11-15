import importlib.util
import unittest
import pandas as pd
from pathlib import Path

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_B_P1_original_syndrome_functional_20261115.py"
spec=importlib.util.spec_from_file_location("alan_syndrome",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestP1Syndrome(unittest.TestCase):
    def test_source_sha_fixed(self):
        self.assertEqual(len(m.SOURCE_SHA),64)
        self.assertEqual(len(m.TRAIT_SHA),64)
        self.assertEqual(m.COMMIT,"1945b7c114ad87000608ce6018d0280001c04ef1")
    def test_no_enforced_complete_cases(self):
        self.assertEqual(m.classify_shape(""),"unknown")
        self.assertEqual(m.classify_shape("brush"),"non_tube")
        self.assertEqual(m.classify_bloom("both"),"both")
    def test_zero_case_study_summary(self):
        x=pd.DataFrame(columns=["study_ID","plant_species","night_loss"])
        d=m.study_summary(x)
        self.assertEqual(d["panels"],0)
        self.assertIsNone(d["study_equal_median"])
    def test_outcome_contrast_no_fictional_3cluster_ci(self):
        data=pd.DataFrame([
          {"study_ID":"a","plant_species":"plant A","night_loss":.8,"trait_profile":"night_tube"},
          {"study_ID":"b","plant_species":"plant B","night_loss":.3,"trait_profile":"day_non_tube"},
          {"study_ID":"c","plant_species":"plant C","night_loss":.2,"trait_profile":"day_non_tube"}
        ])
        c=m.contrast(data,"night_tube","day_non_tube")
        self.assertIsNone(c["paper_bootstrap_95_descriptive"])
        self.assertEqual(c["a"]["study_IDs"],1)
        self.assertIsNotNone(c["median_study_equal_difference"])
    def test_source_paper_loo_not_individual_panel(self):
        src=[]
        for study,value,group in (
          ("N1",.9,"night_tube"),("N2",.8,"night_tube"),("N3",.5,"night_tube"),
          ("D1",.2,"day_non_tube"),("D2",.3,"day_non_tube"),("D3",.4,"day_non_tube"),
        ):
            src.append({"study_ID":study,"plant_species":"plant "+study,
                        "night_loss":value,"trait_profile":group})
        # Extra within-study source panels must not become new paper replicates
        src.append({"study_ID":"N1","plant_species":"plant N1",
                    "night_loss":.95,"trait_profile":"night_tube"})
        c=m.contrast(pd.DataFrame(src),"night_tube","day_non_tube")
        self.assertEqual(c["a"]["study_IDs"],3)
        self.assertEqual(c["b"]["study_IDs"],3)
        self.assertEqual(c["leave_one_paper_n_total"],6)
        self.assertEqual(c["leave_one_paper_n_positive"],6)
        self.assertTrue(c["leave_one_paper_not_proof_independence"])
    def test_original_light_pooled_study_exclusion_label(self):
        self.assertTrue("Macgregor_2019" in SRC.read_text(encoding="utf-8"))
        self.assertTrue("source_removed_panels" in SRC.read_text(encoding="utf-8"))
    def test_strict_source_endpoints(self):
        self.assertIn("seed_per_fruit",m.ENDPOINTS)
        self.assertIn("fruit_proportion",m.ENDPOINTS)
        self.assertNotIn("visit_number",m.ENDPOINTS)

if __name__=="__main__":unittest.main()
