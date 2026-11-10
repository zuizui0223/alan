import unittest
import pandas as pd
from scripts.acquisition.alan_tier09_vs_measured24576_geography_20261110 import score

class MeasuredGeographicSupport(unittest.TestCase):
    def data(self):
        old=pd.DataFrame([{'photo_id':'o1','land_id':'A','year':'2024',
                           'target_id':'2024|1_1|A','strict9':True}])
        new=pd.DataFrame([{'photo_id':f'n{i}','land_id':'A','year':'2024',
                           'cell_id':'2_2','species':('a','b','c')[i%3],
                           'area_km2':12,'viirs_eligible':True} for i in range(6)])
        pairs=pd.DataFrame([{'cohort':'GloBI_moth_record_positive','land_id':'A','year':'2024'}])
        return old,new,pairs

    def test_new_site_not_new_island(self):
        a,b,c=self.data()
        r=score(a,b,c)
        self.assertEqual(r['tier09_supported_islands_without_prior_numeric_ROI9'],0)
        self.assertEqual(r['tier09_qualified_cells_without_prior_numeric_ROI9'],1)
        self.assertEqual(r['tier09_small_islands_le20km2'][0]['tier09_qualified_cells_not_previously_numeric'],1)

    def test_photo_id_collision_fails(self):
        a,b,c=self.data()
        b.loc[0,'photo_id']='o1'
        with self.assertRaisesRegex(ValueError,'original photo ID duplicated'):
            score(a,b,c)

    def test_unknown_is_not_a_success(self):
        a,b,c=self.data()
        a.loc[0,'strict9']=False
        self.assertEqual(score(a,b,c)['tier09_supported_islands_without_prior_numeric_ROI9'],1)

if __name__=='__main__':
    unittest.main()
