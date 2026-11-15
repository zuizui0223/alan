import unittest
import pandas as pd
from scripts.analysis import alan_functional_dependence_20261115 as m

class SinglePaperInterval(unittest.TestCase):
    def test_one_paper_interval_is_unavailable_not_zero_width(self):
        x=pd.DataFrame({'study_ID':['a','a'],'night_loss':[.1,.3]})
        z=m.cluster_summary(x,iterations=20)
        self.assertAlmostEqual(z['median_study_median'],.2)
        self.assertIsNone(z['bootstrap95_median'])
        self.assertIsNone(z['bootstrap95_mean'])

if __name__=='__main__':unittest.main()
