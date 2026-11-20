"""Check block influence, coefficient decomposition and matched ID accounting."""
import importlib.util, pathlib, sys, unittest
import numpy as np
import pandas as pd
S=pathlib.Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(S))
spec=importlib.util.spec_from_file_location("blockdiag",S/"alan_wet_strict_block_influence_20261120.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestBlock(unittest.TestCase):
    def sample(self):
        rng=np.random.default_rng(40)
        t=m.TERMS
        z=pd.DataFrame(rng.normal(size=(90,len(t))),columns=t)
        z.light_difference=np.abs(z.light_difference)+1.
        z["block"]=np.repeat(np.arange(15),6)
        z["candidate_difference"]=-.2*z.light_difference+z.d_built*.1+rng.normal(scale=.1,size=len(z))
        z["reference_difference"]=.1*z.light_difference+rng.normal(scale=.1,size=len(z))
        z["outcome_difference"]=z.candidate_difference-z.reference_difference
        z["bright_site_id"]=np.arange(90)
        z["comparison_site_id"]=np.arange(90)+1000
        return z
    def test_component_decomposition_under_repeated_blocks(self):
        p=self.sample()
        x=m.block_bootstrap_components(p,n_draws=700)
        self.assertEqual(x["status"],"DESCRIPTIVE_COMPONENTS")
        self.assertTrue(x["outcome_identity_verified"])
        self.assertEqual(x["bootstrap_complete_rank_draws"],700)
    def test_loo_all_blocks(self):
        r=m.leave_one_spatial_block_out(self.sample())
        self.assertEqual(r["n_total_blocks"],15)
        self.assertEqual(r["n_fitted"],15)
        self.assertLess(r["max_slope"],0)
    def test_original_pair_identity_accounting(self):
        a=self.sample()
        b=a.iloc[2:75].copy()
        q=m.cohort_relation(a,b)
        self.assertEqual(q["identical_original_pair_id_combinations"],73)
        self.assertEqual(q["current_sites_not_in_original"],0)
    def test_actual_source_cohort_counts_and_pairing(self):
        source,sha=m.read_source()
        self.assertEqual(sha,m.SOURCE_SHA256)
        orig,_,_=m.pair_sites(source,False)
        self.assertEqual((len(orig),orig.block.nunique()),(187,30))
        p,_,_,_=m.stricter_pairs(source,.5)
        q,_,_,_=m.stricter_pairs(source,.25)
        self.assertEqual((len(p),p.block.nunique()),(133,26))
        self.assertEqual((len(q),q.block.nunique()),(128,25))

if __name__=="__main__":
    unittest.main()
