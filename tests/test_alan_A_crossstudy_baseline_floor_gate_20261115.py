import json,unittest,copy
from scripts.analysis import alan_A_crossstudy_baseline_floor_gate_20261115 as m
from scripts.analysis import alan_Grognuz2026_block_floor_signs_20261115 as grog
class TestTwoStudyGate(unittest.TestCase):
    def test_author_paper_ids_and_experimental_k(self):
        reg=json.loads(m.SOURCE_REGISTER.read_text(encoding="utf-8"))
        eligible={x["id"] for x in reg["entries"]
          if x["light_design"].startswith("manipulated")
          and x["a_primary_stage"]!="exclude_A_no_reproduction"}
        self.assertEqual(len(eligible),9)
    def test_grognuz_block_mean_and_zero(self):
        d=grog.analyze(grog.load())
        self.assertEqual(d["total_blocks"],28)
        self.assertEqual(d["Trifolium_pratense"]["dark_zero_blocks"],["T03"])
        self.assertEqual(d["Centaurea_jacea"]["dark_zero_blocks"],["C01"])
    def test_crossstudy_never_claims_k39(self):
        self.assertIn("true_independent_lighting_sites_not_11_plus_28",__import__("pathlib").Path(m.__file__).read_text(encoding="utf-8"))
    def test_covariance_identity(self):
        import numpy as np
        a=np.asarray([0.,1.,2.,3.]);b=np.asarray([2.,3.,1.,5.])
        self.assertAlmostEqual(float(np.cov(a,b-a,bias=True)[0,1]),
            float(np.cov(a,b,bias=True)[0,1]-np.var(a)))
if __name__=="__main__":unittest.main()
