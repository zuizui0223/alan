"""Synthetic source-only math checks, not new ecological observations."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'/'analysis'))
import alan_13year_moth_gecko_repeat_site_FE_20261112 as m

def generated():
    rand=np.random.default_rng(456)
    out=[]
    for island in range(15):
        for year in (2019,2020,2021,2022):
            for cell in range(4):
                loglight=float(.3+.04*(year-2019)+.3*cell+.05*rand.normal())
                gecko=float(rand.beta(1.7,1.7))
                log_insect=float(rand.uniform(1,5))
                log_reptile=float(rand.uniform(1,3))
                noise=float(rand.normal(0,.3))
                result=(1.25*loglight+2.8*gecko-2.3*loglight*gecko
                         +.6*log_insect-.2*log_reptile+noise
                         +.3*island+.1*(year-2019))
                out.append(dict(island=island,year=year,
                    site=f'isl{island}|grid{cell}',
                    islandyear=f'isl{island}|{year}',
                    loglight=loglight,gecko_share=gecko,
                    cross=loglight*gecko,
                    log_insect=log_insect,log_reptile=log_reptile,
                    moth_share=result,
                    moth_family_records=5,gecko_source_photo_count=2))
    return pd.DataFrame(out)

class SourceFixedEffects(unittest.TestCase):
    def test_sparse_projection_equals_explicit_weighted_dummy_WLS(self):
        d=generated()
        fit=m.fit(d)
        self.assertTrue(fit['estimable'])
        self.assertEqual(fit['source_connected_FE_components'],15)
        self.assertEqual(fit['fixed_effects'],{'site':60,'islandyear':60})
        w=1/d.groupby('island').island.transform('size').to_numpy(float)
        reference=smf.wls('moth_share ~ loglight + gecko_share + cross + log_insect + log_reptile + C(site) + C(islandyear)',data=d,weights=w).fit()
        for term in ('loglight','gecko_share','cross','log_insect','log_reptile'):
            self.assertAlmostEqual(fit['coefficients'][term]['beta'],
                                   reference.params[term],places=8)
        self.assertEqual(fit['effective_residual_df'],len(d)-(60+60-15)-5)

    def test_exact_repeated_site_not_other_island_or_single_year(self):
        d=pd.DataFrame([
          {'site':'A|5_5','year':2019,'loglight':.3},
          {'site':'A|5_5','year':2020,'loglight':.35},
          {'site':'A|6_6','year':2019,'loglight':.1},
          {'site':'B|5_5','year':2020,'loglight':.8},
          {'site':'B|5_5','year':2021,'loglight':.801},
        ])
        z=m.repeated_sites(d,.02)
        self.assertEqual(len(z),2)
        self.assertEqual(z.site.nunique(),1)
        self.assertEqual(z.site.iloc[0],'A|5_5')

    def test_sparse_effects_finite_and_positive_source_only(self):
        d=generated()
        for w in ('cell_equal','island_equal','islandyear_equal'):
            result=m.fit(d,weight=w)
            self.assertTrue(result['estimable'])
            self.assertTrue(np.isfinite(result['coefficients']['cross']['beta']))
            self.assertEqual(result['independent_islands'],15)

    def test_unestimable_in_synthetic_too_few_islands(self):
        d=generated().loc[lambda x:x.island.lt(12)].copy()
        self.assertFalse(m.fit(d)['estimable'])

if __name__=='__main__':
    unittest.main()
