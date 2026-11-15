import importlib.util
import math
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_A_shared_control_null_20261115.py"
SPEC=importlib.util.spec_from_file_location("alan_a_null",SRC)
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

class TestSharedControlNull(unittest.TestCase):
    def test_negative_covariance_under_shared_open(self):
        for day in (0,.1,.2,.5):
            value=M.delta_shared_covariance(day,.6,60)
            self.assertLessEqual(value,0)
        self.assertEqual(M.delta_shared_covariance(0,.6,60),0)
    def test_known_covariance_number(self):
        # Var(O)=0.6*0.4/60, D=0.3
        self.assertAlmostEqual(M.delta_shared_covariance(.3,.6,60),-.3*(.6*.4/60)/(.6**3))
    def test_no_effect_still_negative_shared(self):
        s=M.simulate_one(60,20261175,replicates=10000)
        self.assertLess(s["shared_dark_open"]["mean_beta"],-0.03)
        self.assertLess(abs(s["independent_moderator_open"]["mean_beta"]),.05)
        self.assertEqual(s["true_ALAN_log_RR"],0)
    def test_reproducibility(self):
        a=M.simulate_one(60,20261175,replicates=200)
        b=M.simulate_one(60,20261175,replicates=200)
        self.assertEqual(a,b)

if __name__=="__main__":
    unittest.main()
