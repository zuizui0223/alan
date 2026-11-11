"""Synthetic source safeguards and mathematical checks; not ecological observations."""
import gzip
import hashlib
import io
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts' / 'analysis'))
import alan_globi_fiveyear_temporal as temporal
import alan_globi_photo_month_source as season

def make_sources(seed=101):
    rng=np.random.default_rng(seed)
    out=[]
    for year in temporal.YEARS:
        for region in range(15):
            for land in ('island','mainland'):
                for k in range(6):
                    n=int(rng.integers(3,9))
                    m=int(rng.binomial(n,.3))
                    b=int(rng.binomial(n,.4))
                    out.append({
                       'target_id':f'{year}|{region*100+k}_1400|{land}',
                       'land_class':land,'region5':f'reg{region}','year':year,
                       'n_species':n,'photo_attempts':n+1,
                       'viirs_median':float(rng.uniform(.1,8)),
                       'pop_density':float(rng.uniform(.1,400)),
                       'bio5':float(rng.normal(21,3)),
                       'rain':float(rng.uniform(100,3000)),
                       'elevation':float(rng.uniform(0,1000)),
                       'all_moth':m,'all_butterfly':b,
                       'non_inat_moth':max(0,m-1),
                       'non_inat_butterfly':max(0,b-1)})
    return pd.DataFrame(out)

class FiveYearTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.csv=self.root/'PRIVATE_fixture.csv'
        make_sources().to_csv(self.csv,index=False)
        self.frame=temporal.select(temporal.load(self.csv))

    def test_full_years_and_repeated_site_identity(self):
        self.assertEqual(len(self.frame),900)
        self.assertEqual(self.frame.site.nunique(),180)
        self.assertEqual(len(temporal.repeated(self.frame)),900)

    def test_exact_weighted_FWL_equals_explicit_WLS(self):
        x=self.frame.copy()
        model=temporal.fit(x,temporal.basic_terms(),'share_all_moth')
        weights=temporal.balanced_weights(x)
        reference=smf.wls('share_all_moth ~ island + loglight + light_island + logpop + logphotos + logspecies + bio5 + lograin + elevation_km + C(region_year)',data=x,weights=weights).fit()
        self.assertTrue(model['estimable'])
        for term in ('light_island','loglight','island','logphotos'):
            self.assertAlmostEqual(model['coefficients'][term]['beta'],
                                   reference.params[term],places=8)

    def test_moth_minus_butterfly_year_trend_algebra(self):
        terms=temporal.trend_terms()
        m=temporal.fit(self.frame,terms,'share_all_moth')
        b=temporal.fit(self.frame,terms,'share_all_butterfly')
        diff=temporal.fit(self.frame,terms,'share_moth_minus_butterfly')
        for term in ('light_island','year_island_light','year_light'):
            self.assertAlmostEqual(diff['coefficients'][term]['beta'],
                m['coefficients'][term]['beta']-b['coefficients'][term]['beta'],places=8)

    def test_year_heterogeneity_and_repeated_site_FE(self):
        fit=temporal.fit(self.frame,temporal.dynamic_terms(),'share_all_moth')
        test=temporal.year_heterogeneity(fit)
        self.assertTrue(test['estimable'])
        self.assertTrue(0<=test['p_exploratory']<=1)
        terms={'loglight':'loglight','light_island':'light_island',
            'logpop':'logpop','logphotos':'logphotos','logspecies':'logspecies'}
        result=temporal.fit(temporal.repeated(self.frame),terms,'share_all_moth',
                            fe=('site','region_year'),site_fe_df=True)
        self.assertTrue(result['estimable'])

    def test_invalid_species_count_or_calendar_year_rejected(self):
        d=make_sources();d.loc[0,'all_moth']=d.loc[0,'n_species']+1
        d.to_csv(self.csv,index=False)
        with self.assertRaisesRegex(ValueError,'exceeds'):
            temporal.load(self.csv)
        d.loc[0,'all_moth']=0
        d.loc[0,'year']=2022
        d.to_csv(self.csv,index=False)
        with self.assertRaisesRegex(ValueError,'five-year window'):
            temporal.load(self.csv)

    def test_photo_season_target_mismatch_rejected(self):
        d=make_sources()
        season=pd.DataFrame({'target_id':d.target_id.iloc[:-1],
            'valid_dates':1,'sin_month':0.,'cos_month':1.,
            'south_sin':0.,'south_cos':0.})
        seasonal_file=self.root/'PRIVATE_season.csv'
        season.to_csv(seasonal_file,index=False)
        with self.assertRaisesRegex(ValueError,'not aligned'):
            temporal.run(self.csv,seasonal_file)

class PhotoMonthTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.zip=self.root/'source.zip'
        self.sites=self.root/'PRIVATE_sites.csv'
        records=[]
        targets=[]
        for y in temporal.YEARS:
            targets.append({'target_id':f'{y}|3601_1801|LAND',
                            'year':y,'photo_attempts':2})
            for j,month in enumerate((1,7)):
                records.append({'photo_id':f'{y}_{j}',
                    'observed_on':f'{y}-{month:02d}-15',
                    'source_subarm':'historical_main',
                    'latitude':.06,'longitude':.06})
        raw=gzip.compress(pd.DataFrame(records).to_csv(index=False).encode())
        self.source_sha=hashlib.sha256(raw).hexdigest()
        with zipfile.ZipFile(self.zip,'w') as z:
            z.writestr(season.PHOTO_MEMBER,raw)
        pd.DataFrame(targets).to_csv(self.sites,index=False)

    def test_all_months_match_original_photo_opportunities(self):
        private=self.root/'PRIVATE_season.csv.gz'
        public=self.root/'public_receipt.json'
        r=season.run(self.zip,self.sites,private,public,
                     expected_sha=self.source_sha,synthetic_fixture=True)
        self.assertEqual(r['matched_original_photo_attempts'],10)
        self.assertEqual(r['site_years_exact_photo_attempts_reconciled'],5)
        self.assertNotIn('longitude',public.read_text())
        self.assertEqual(len(pd.read_csv(private)),5)

    def test_source_count_mismatch_refused(self):
        d=pd.read_csv(self.sites)
        d.loc[0,'photo_attempts']=3
        d.to_csv(self.sites,index=False)
        with self.assertRaisesRegex(ValueError,'disagree'):
            season.run(self.zip,self.sites,self.root/'PRIVATE_season.csv.gz',
                       self.root/'receipt.json',
                       expected_sha=self.source_sha,synthetic_fixture=True)

    def test_changed_canonical_SHA_refused(self):
        with self.assertRaisesRegex(ValueError,'SHA mismatch'):
            season.run(self.zip,self.sites,self.root/'PRIVATE_season.csv.gz',
                       self.root/'receipt.json',
                       expected_sha='0'*64,synthetic_fixture=True)

if __name__=='__main__':
    unittest.main()
