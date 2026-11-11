"""Synthetic regression invariance checks. No fictional ecological data published."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_globi_crossland_light_robust_20261111 as model

def fixture():
    rg=np.random.default_rng(101)
    rows=[]
    for region in range(18):
        for year in (2019,2023):
            for land in ('island','mainland'):
                for k in range(12):
                    light=rg.uniform(.1,4)
                    n=int(rg.integers(3,13))
                    moth=int(rg.binomial(n,.25+.06*light/(1+light)+.04*(land=='island')))
                    butterfly=int(rg.binomial(n,.4+.02*light/(1+light)))
                    row={'target_id':f'{year}_{region}_{land}_{k}','land_class':land,
                         'n_species':n,'viirs_median':light,'photo_attempts':n+int(rg.integers(0,6)),
                         'bio5':rg.normal(20,3),'rain':rg.uniform(500,2000),
                         'elevation':rg.uniform(0,2000),'pop_density':rg.uniform(.1,500),
                         'year':year,'region5':f'block_{region}'}
                    for scope in model.SCOPES:
                        row[f'{scope}_moth']=moth
                        row[f'{scope}_butterfly']=butterfly
                    rows.append(row)
    return pd.DataFrame(rows)

class TestGloBICrossland(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.file=Path(self.temp.name)/'PRIVATE_fake.csv'
        self.data=fixture()
        self.data.to_csv(self.file,index=False)
        self.frame=model.cleaned(self.file)

    def test_contrast_is_exact_difference_of_two_coefficients(self):
        matched=model.prepare(self.frame,3,2)
        for balanced in (False,True):
            moth=model.fit(matched,'all','moth',balance=balanced)
            butterfly=model.fit(matched,'all','butterfly',balance=balanced)
            contrast=model.fit(matched,'all','moth_minus_butterfly',balance=balanced)
            self.assertTrue(all(x['estimable'] for x in (moth,butterfly,contrast)))
            self.assertAlmostEqual(contrast['beta_island_minus_mainland_slope_pp_per_log1p'],
                moth['beta_island_minus_mainland_slope_pp_per_log1p'] -
                butterfly['beta_island_minus_mainland_slope_pp_per_log1p'],places=8)

    def test_no_positive_globi_record_yields_zero_evidence_not_absence(self):
        frame=self.frame.copy()
        frame['non_inat_moth']=0
        fit=model.fit(model.prepare(frame,3,2),'non_inat','moth',balance=True)
        self.assertAlmostEqual(fit['beta_island_minus_mainland_slope_pp_per_log1p'],0,places=8)
        self.assertTrue(fit['no_true_pollinator_absence_inferred'])

    def test_disjoint_light_ranges_fail_overlap_screen(self):
        frame=self.frame.copy()
        frame.loc[frame.land_class.eq('island'),'viirs_median']=50
        frame.loc[frame.land_class.eq('mainland'),'viirs_median']=.02
        frame.to_csv(self.file,index=False)
        clean=model.cleaned(self.file)
        self.assertEqual(len(model.prepare(clean,3,2,light_overlap=True)),0)

    def test_bad_species_count_fails_closed(self):
        frame=self.data.copy()
        frame.loc[0,'all_moth']=frame.loc[0,'n_species']+1
        frame.to_csv(self.file,index=False)
        with self.assertRaisesRegex(ValueError,'count exceeds'):
            model.cleaned(self.file)

    def test_support_and_fits(self):
        result=model.run(self.file)
        self.assertEqual(len(result['models']),60)
        left=result['source_support']['3photos_2species']['class_counts']['island']['n']
        right=result['source_support']['3photos_2species']['class_counts']['mainland']['n']
        self.assertEqual(left,right)

if __name__=='__main__':
    unittest.main()
