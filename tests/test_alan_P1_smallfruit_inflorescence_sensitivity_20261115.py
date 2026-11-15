import importlib.util
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_P1_smallfruit_inflorescence_sensitivity_20261115.py"
spec=importlib.util.spec_from_file_location("alan_p1_b1",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def row(taxon,year,arm,success,failed):
    return {"Plant_sp":taxon,"Year":str(year),"Treatment":arm,
      "Fruit_set_success":str(success),"Fruit_set_failure":str(failed)}

class TestSmallFruitSensitivity(unittest.TestCase):
    def test_original_data_inference_unit(self):
        self.assertEqual(m.MAP["BC"],"Ribes nigrum")
        self.assertEqual(m.MAP["RB"],"Rubus idaeus")
        self.assertEqual(m.MAP["RD"],"Ribes rubrum")
    def test_n_zero_open_no_invented_values(self):
        day=[(2,8),(3,7)]
        no_open=[(0,10),(0,10)]
        x=m.diagnostics(day,no_open,2026,boot_n=20)
        self.assertIsNone(x["V"])
    def test_small_2020_arms_no_fake_CI(self):
        day=[(1,9),(3,7)]
        op=[(5,5),(5,5)]
        x=m.diagnostics(day,op,42,boot_n=10)
        self.assertAlmostEqual(x["V"],.6)
        self.assertIsNone(x["CI_descriptive"])
        self.assertEqual(x["ci_status"],"SUPPRESSED_LT3_INFLORESCENCES_PER_ARM")
        self.assertEqual(x["leave_one_signs"]["positive"],4)
    def test_leave_one_can_detect_reverse_sign(self):
        day=[(8,2),(0,10),(0,10)]
        op=[(2,8),(2,8),(2,8)]
        x=m.diagnostics(day,op,42,boot_n=100)
        self.assertLess(x["leave_one_inflorescence_min"],0)
        self.assertGreater(x["leave_one_inflorescence_max"],0)
        self.assertEqual(x["true_plant_cluster_CI"],None)
    def test_reproducible_bootstrap_fixed_seed(self):
        d=[(1,7),(2,8),(2,7)]
        o=[(4,6),(5,5),(6,4)]
        a=m.diagnostics(d,o,42,boot_n=120)
        b=m.diagnostics(d,o,42,boot_n=120)
        self.assertEqual(a,b)
        self.assertIsNotNone(a["CI_descriptive"])
    def test_full_four_arm_from_synthetic_plant_data(self):
        rows=[
          row("BC",2024,"Closed",1,9),row("BC",2024,"Day",2,8),
          row("BC",2024,"Night",4,6),row("BC",2024,"Day + Night",5,5),
        ]
        out=m.analyze(rows,boot_n=100)
        r=out["results"][0]
        self.assertAlmostEqual(r["night_access_marginal_V"]["V"],.6)
        self.assertIsNone(r["night_access_marginal_V"]["CI_descriptive"])
        self.assertAlmostEqual(r["closed_over_open"],.2)
        self.assertFalse(out["moth_pollinator_identification"])
        self.assertTrue(out["no_ALAN_manipulated"])

if __name__ == "__main__":
    unittest.main()
