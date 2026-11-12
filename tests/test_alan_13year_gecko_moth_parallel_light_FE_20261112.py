"""Synthetic independent WLS math checks for parallel gecko/moth photo-ratio models."""
import sys,unittest
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_13year_gecko_moth_parallel_light_FE_20261112 as m

def synthetic():
    rng=np.random.default_rng(919)
    z=[]
    for i in range(15):
        for y in (2019,2020,2021,2022):
            for j in range(4):
                L=.4+.07*(y-2019)+.25*j+rng.normal(0,.06)
                s=rng.beta(2,2)
                pp=rng.uniform(1,4);qq=rng.uniform(1,3)
                z.append({'island':i,'year':y,
                 'site':f'i{i}|cell{j}','islandyear':f'i{i}|{y}',
                 'loglight':L,'gecko_share':s,
                 'moth_share':80-2*L+2*s+rng.normal(0,2),
                 'log_insect':pp,'log_reptile':qq,
                 'latitude':20+.3*i+.01*j})
    return pd.DataFrame(z)

class MarginalSourceModels(unittest.TestCase):
    def test_spatial_against_explicit_islandyear_WLS(self):
        z=synthetic()
        for guild in ('gecko','moth'):
            model=m.fit(z,guild)
            self.assertTrue(model['estimable'])
            z['outcome']=(100*z.gecko_share) if guild=='gecko' else z.moth_share
            w=1/z.groupby('island').island.transform('size').to_numpy(float)
            ref=smf.wls('outcome ~ loglight + log_insect + log_reptile + latitude + C(islandyear)',data=z,weights=w).fit()
            self.assertAlmostEqual(model['slope_light'],ref.params['loglight'],places=8)

    def test_same_site_temporal_against_explicit_two_way_WLS(self):
        z=synthetic()
        for guild in ('gecko','moth'):
            model=m.fit(z,guild,True)
            self.assertTrue(model['estimable'])
            self.assertEqual(model['FE_connected_components'],15)
            z['outcome']=(100*z.gecko_share) if guild=='gecko' else z.moth_share
            w=1/z.groupby('island').island.transform('size').to_numpy(float)
            ref=smf.wls('outcome ~ loglight + log_insect + log_reptile + C(site) + C(islandyear)',data=z,weights=w).fit()
            self.assertAlmostEqual(model['slope_light'],ref.params['loglight'],places=8)

    def test_photo_observation_design_requires_geographic_clusters(self):
        z=synthetic().loc[lambda x:x.island.lt(14)]
        self.assertFalse(m.fit(z)['estimable'])

if __name__=='__main__':
    unittest.main()
