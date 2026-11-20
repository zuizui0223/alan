"""Secondary per-mesocosm period means, keeping final zero arms visible."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/bennie2018_trophic'
manifest=[]
for f in json.loads((CACHE/'files.json').read_text())['_embedded']['stash:files']:
    p=CACHE/f['path'];digest=hashlib.new(f['digestType'].replace('-',''),p.read_bytes()).hexdigest()
    assert digest==f['digest']
    manifest.append(dict(file=f['path'],algorithm=f['digestType'],digest=digest,dryad_match=True))
a=pd.read_csv(CACHE/'slugs2013.csv');b=pd.read_csv(CACHE/'slugs2014.csv')
keys=['Plot','T.treatment','L.treatment']
months13=['June','July','August','September','October','November']
months14=['January','February','March','April','May','June','July','August','September']
assert len(a)==len(b)==36 and a.Plot.is_unique and b.Plot.is_unique
m=a[keys+months13].merge(b[keys+months14],on=keys,validate='one_to_one',suffixes=('_2013','_2014'))
assert len(m)==36
values=m.drop(columns=keys)
assert values.shape==(36,15) and not values.isna().any().any()
m['period_mean']=values.mean(axis=1)
m['months']=15
m['final_2014_September']=m.September_2014
m.to_csv(OUT/'slug_original_unit_period_means.csv',index=False)
# Physical allocation matches existing Penryn54; no new independent programme.
cross=pd.read_csv(ROOT/'database/alan_experiments/bennie_sources/anic2022_same_mesocosm_crosswalk.csv')
check=m[keys].merge(cross,on=None,left_on='Plot',right_on='plot',validate='one_to_one')
assert len(check)==36
assert (check['T.treatment'].map({'Grazer':'T2','Predator':'T3'})==check.trophic_2013).all()
assert (check['L.treatment'].replace({'Control':'Dark'})==check.light_2013).all()
check[keys+['trophic_2013','light_2013','light_2016','trophic_2016']].to_csv(OUT/'programme_crosswalk.csv',index=False)
# Also confirm repeated 2013 aphid endpoint is already represented, not new.
aph=pd.read_csv(CACHE/'PeaAphid2013.csv')[['Plot','12/09/2013']]
old=pd.read_csv(ROOT/'database/alan_experiments/bennie_sources/final_original_units.csv')
old=old[old.endpoint=='aphid_count_3min']
j=aph.merge(old[['plot','value']],left_on='Plot',right_on='plot',validate='one_to_one')
assert len(j)==54 and np.allclose(j['12/09/2013'],j.value,equal_nan=True)
aph_all=pd.read_csv(CACHE/'PeaAphid2013.csv')
aph_long=aph_all.melt(id_vars=['Plot','T treatment','L treatment'],var_name='date',value_name='new_value')
aph_long['date']=pd.to_datetime(aph_long.date,dayfirst=True).dt.strftime('%Y-%m-%d')
old_long=pd.read_csv(ROOT/'database/alan_experiments/bennie_sources/repeated_original_units.csv')
old_long=old_long[old_long.endpoint=='aphid_count_3min']
comparison=aph_long.merge(old_long[['plot','date','value']],left_on=['Plot','date'],right_on=['plot','date'],validate='one_to_one')
assert len(comparison)==378
differences=comparison[~np.isclose(comparison.new_value,comparison.value,equal_nan=True)]
assert len(differences)==2
differences[['Plot','date','new_value','value']].rename(columns={'new_value':'Bennie2018_value','value':'Bennie2015_value'}).to_csv(OUT/'aphid_2013_source_discrepancies.csv',index=False)
cells=[(g,l) for g in ['Grazer','Predator'] for l in ['Control','Amber','White']]
arms=[]
for g,l in cells:
    z=m[(m['T.treatment']==g)&(m['L.treatment']==l)]
    assert len(z)==6
    arms.append(dict(trophic=g,light=l,n_mesocosms=6,months=15,
                     mean=z.period_mean.mean(),sd=z.period_mean.std(ddof=1),
                     final_mean=z.final_2014_September.mean(),final_sd=z.final_2014_September.std(ddof=1)))
arms=pd.DataFrame(arms);arms.to_csv(OUT/'arm_summaries.csv',index=False)
mu=arms['mean'].to_numpy();var=arms.sd.to_numpy()**2/6
J=np.zeros((4,6));theta=[];ids=[]
for i,(g,l) in enumerate([(g,l) for g in ['Grazer','Predator'] for l in ['Amber','White']]):
    c=cells.index((g,'Control'));k=cells.index((g,l))
    J[i,c]=-1/mu[c];J[i,k]=1/mu[k]
    theta.append(np.log(mu[k]/mu[c]));ids.append(f'Bennie2018_slug_period_{g}_{l}')
theta=np.array(theta);V=J@np.diag(var)@J.T
A=np.array([[-1.,0.,1.,0.],[0.,-1.,0.,1.]])
def save(stem,ids,est,cov):
    half=t.ppf(.975,5)*np.sqrt(np.diag(cov))
    rows=pd.DataFrame(dict(effect_id=ids,lnRR=est,sampling_variance=np.diag(cov),ratio=np.exp(est),
                           approximate_t5_lower=np.exp(est-half),approximate_t5_upper=np.exp(est+half),
                           main_catalogue=False,programme_id='Penryn2012_2016',estimand='secondary_equal_month_period_mean'))
    rows.to_csv(OUT/f'{stem}.csv',index=False)
    pd.DataFrame([dict(effect_id_i=ids[i],effect_id_j=ids[j],covariance=cov[i,j]) for i in range(len(ids)) for j in range(len(ids))]).to_csv(OUT/f'{stem}_covariance.csv',index=False)
    assert np.linalg.eigvalsh(cov).min()>-1e-12
    return rows
save('slug_period_effects',ids,theta,V)
out=save('slug_predator_contrasts',['Bennie2018_slug_period_predator_Amber','Bennie2018_slug_period_predator_White'],A@theta,A@V@A.T)
assert arms.loc[arms.light=='Control','final_mean'].eq(0).all()
# Final endpoint remains untransformed; never replace zero dark means.
final=[]
for g in ['Grazer','Predator']:
    dark=arms[(arms.trophic==g)&(arms.light=='Control')].iloc[0]
    for l in ['Amber','White']:
        lit=arms[(arms.trophic==g)&(arms.light==l)].iloc[0]
        final.append(dict(trophic=g,light=l,mean_lit=lit.final_mean,mean_dark=dark.final_mean,
                          mean_difference=lit.final_mean-dark.final_mean,
                          mean_difference_variance=(lit.final_sd**2+dark.final_sd**2)/6,
                          lnRR_status='undefined_zero_dark_mean',n_lit=6,n_dark=6))
pd.DataFrame(final).to_csv(OUT/'final_2014_September_zero_arm_results.csv',index=False)
audit=dict(files=manifest,mirror='https://zenodo.org/records/4988926',paper_doi='10.1111/1365-2664.13240',
           slug_mesocosms=36,months_per_unit=15,independent_units_per_cell=6,
           period='June2013-November2013 and January2014-September2014; December2013 not supplied',
           zero_final_dark_groups=2,finite_primary_final_lnRR=0,secondary_period_lnRR=4,
           dependent_predator_contrasts=2,matched_existing_programme_plots=36,
           repeated_2013_aphid_endpoint_rows_matched=54,new_independent_programmes=0,
           intermediate_2013_aphid_value_discrepancies=2,compared_2013_aphid_plot_dates=378,
           covariance_structure='independent randomized mesocosm arms; shared dark arms; not paired by plot letter',
           main_catalogue=False,interpretation='Counts, not attack rate; equal-month period mean is not author longitudinal GLMM.')
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(out[['effect_id','ratio','approximate_t5_lower','approximate_t5_upper']].to_string(index=False))
