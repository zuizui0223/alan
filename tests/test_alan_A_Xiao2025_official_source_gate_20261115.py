import json
from pathlib import Path
import unittest

REG=json.loads((Path(__file__).resolve().parents[1]/"config/alan_A_reproductive_screen_20261115.json").read_text(encoding="utf-8"))
X=next(x for x in REG["entries"] if x["id"]=="xiao2025")

class TestXiao2025OriginalReproductiveSource(unittest.TestCase):
 def test_real_40plot_treatment_factorial_not_moth_exclosure(self):
    a=X["source_audit"]
    self.assertEqual(a["experimental_plot_count"],40)
    self.assertEqual(a["plot_area_m2"],4)
    self.assertEqual(a["factorial"],"ALAN/night dark x animal exclusion/open")
    self.assertIn("not authenticated",a["animal_exclusion_mechanism"])
 def test_actual_official_supplement_bytes_and_missing_tables(self):
    a=X["source_audit"]
    self.assertEqual(a["official_supplement_http"],200)
    self.assertEqual(a["supplementary_file_size_bytes"],315350)
    self.assertEqual(a["supplementary_file_sha256"],"127950cecd9ac78c0927b89df06b307e92a7a9445b0c054aa40a2e9d272e8a9f")
    self.assertEqual(a["supplementary_docx_table_count"],0)
    self.assertEqual(a["supplementary_embedded_image_count"],3)
    self.assertFalse(a["source_reproduction_point_estimate_and_se_available"])
 def test_reproduction_is_not_initial_flower_success(self):
    self.assertIn("secondary_fruiting_inflorescence_count",X["a_primary_stage"])
    self.assertNotIn("successful_fruit_per_initial_flower",X["outcome"])
    self.assertFalse(X["verified_joint_effect_moderator_se"])
    self.assertFalse(X["source_audit"]["night_insect_pollination_dependency_V_available"])
 def test_candidate_count_unchanged(self):
    self.assertEqual(REG["design_counts_among_reproductive"],{"manipulated":9,"observational":3,"mixed_unresolved":0})
if __name__=="__main__":unittest.main()
