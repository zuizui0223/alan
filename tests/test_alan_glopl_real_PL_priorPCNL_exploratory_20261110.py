"""Synthetic support-gates only: PL effect requires independent original experiment evidence."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_glopl_real_PL_priorPCNL_exploratory_20261110 as m


class TestGloPLAssociationGuards(unittest.TestCase):
    def test_original_PL_cannot_be_stand_in_for_moth_pollination(self):
        self.assertEqual(len(m.METH_COLS),4)
        self.assertEqual(m.EXPOSURE,"pcnl_prior_mean_log1p")
        self.assertIn("PL",m.MODEL_SCHEMA)
    def test_insufficient_publication_support_not_estimated(self):
        x=pd.DataFrame({"study_key":["a","a","b"],m.EXPOSURE:[0.,1.,2.],
            "publication_equal_weight":[.5,.5,1.],
            "lat":[0.,0.,0.],"lon":[0.,1.,2.],"year_exp":[2010,2010,2010],
            "PL_log_ratio":[.4,.2,.8]})
        f=m.fit(x,"PL_log_ratio ~ pcnl_prior_mean_log1p","study_key")
        self.assertEqual(f["status"],"insufficient_source_support")
    def test_source_prior_year_calendrical_source_is_bound_to_experiment(self):
        self.assertEqual(m.SOURCE_SCHEMA,"ALAN_real_GloPL_PCNLPrior1992_2014_source_join_v1")


if __name__=="__main__":unittest.main()
