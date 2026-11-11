"""Synthetic source geometry and WLS controls. NOT ecological evidence."""
import hashlib,sys,tempfile,unittest
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_globi_crossguild_photo_light_replay_20261111 as m

def example():
    rand=np.random.default_rng(42)
    rows=[]
    for year in (2019,2024):
        for region in range(17):
            for land in ('island','mainland'):
                for k in range(12):
                    n=int(rand.integers(3,12))
                    one={'target_id':f'{year}|{region}_{k}|{land}',
                      'land_id':f'{region}_{land}','year':year,'region5':f'r{region}',
                      'land_class':land,'viirs':float(rand.uniform(.01,25)),
                      'population_density':float(rand.uniform(.1,300)),
                      'photo_attempts':n+1,'n_species':n,
                      'temperature':float(rand.normal(19,3)),
                      'precipitation':float(rand.uniform(200,1800)),
                      'elevation':float(rand.uniform(1,1500))}
                    for scope in m.SCOPES:
                        for group in m.GROUPS:
                            for relation in m.RELATIONS:
                                one[f'count_{scope}_{group}_{relation}']=(
                                  int(rand.binomial(n,.28)) if scope=='all' and group=='birds'
                                  else 0)
                    rows.append(one)
    return pd.DataFrame(rows)

class Tests(unittest.TestCase):
    def test_source_sha_and_missing_as_unrecorded_not_biological_absence(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'private.csv.gz'
            a=example()
            a.to_csv(path,index=False,compression={'method':'gzip','mtime':0})
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            d=m.run(path,expected_sha=digest)
            birds=next(x for x in d['models'] if x['source_scope']=='all' and x['animal_guild']=='birds' and x['GloBI_relation']=='flower_either')
            geckos=next(x for x in d['models'] if x['animal_guild']=='geckos')
            self.assertTrue(birds['test']['estimable'])
            self.assertFalse(geckos['test']['estimable'])
            self.assertTrue(d['GloBI_global_visitor_positive_not_local_pollination'])
            with self.assertRaisesRegex(ValueError,'SHA256 mismatch'):
                m.load(path,expected_sha='0'*64)

    def test_fixed_effect_regression_equals_independent_WLS(self):
        data=example()
        data['island']=data.land_class.eq('island').astype(int)
        data['region_year']=data.year.astype(str)+'|'+data.region5
        response='count_all_birds_flower_either'
        fitted=m.fit(data,response)
        self.assertTrue(fitted['estimable'])
        d=data.copy()
        d['response']=100*d[response]/d.n_species
        d['L']=np.log1p(d.viirs)
        d['logpop']=np.log1p(d.population_density)
        d['logphotos']=np.log1p(d.photo_attempts)
        d['logspecies']=np.log1p(d.n_species)
        d['lograin']=np.log1p(d.precipitation)
        d['elevation_km']=d.elevation/1000
        n=d.groupby(['region_year','land_class']).size()
        ix=pd.MultiIndex.from_frame(d[['region_year','land_class']])
        w=n.groupby(level=0).min().reindex(d.region_year).to_numpy(float)/n.reindex(ix).to_numpy(float)
        reference=smf.wls('response ~ L*island + logpop + logphotos + logspecies + temperature + lograin + elevation_km + C(region_year)',d,weights=w).fit()
        self.assertAlmostEqual(fitted['beta_island_minus_mainland_pp_per_log1pVIIRS'],reference.params['L:island'],places=8)

    def test_one_sided_year_region_not_used(self):
        d=example()
        d=d.loc[~((d.year==2024)&(d.region5=='r0')&(d.land_class=='island'))].copy()
        d['island']=d.land_class.eq('island').astype(int)
        d['region_year']=d.year.astype(str)+'|'+d.region5
        result=m.fit(d,'count_all_birds_flower_either')
        self.assertTrue(result['estimable'])
        # Year r0/2024 mainland-only is excluded. 17*2*24 - 12 retained
        self.assertEqual(result['source_site_years'],(17*2-1)*24)

    def test_duplicates_and_positive_count_excess_fail(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'private.csv.gz'
            d=example()
            d.loc[0,'count_all_bats_flower_either']=d.loc[0,'n_species']+1
            d.to_csv(path,index=False,compression='gzip')
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError,'exceed'):
                m.load(path,digest)

if __name__=='__main__':unittest.main()
