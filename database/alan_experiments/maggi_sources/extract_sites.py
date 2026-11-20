"""Site-level extraction; pigment composition is not biomass or grazing rate."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/maggi2020'
manifest = []
for f in json.loads((CACHE/'files.json').read_text())['_embedded']['stash:files']:
    p = CACHE/f['path']
    actual = hashlib.new(f['digestType'].replace('-',''),p.read_bytes()).hexdigest()
    assert actual == f['digest']
    manifest.append(dict(file=f['path'],algorithm=f['digestType'],digest=actual,dryad_match=True))
snails = pd.read_excel(CACHE/'data_Littorinid_snails.xlsx')
assert len(snails)==216 and not snails.duplicated(['site','plot','time']).any()
assert snails.groupby(['site','plot']).time.nunique().eq(9).all()
snail_sites = snails.groupby(['light','site','herb','time'],as_index=False).agg(mean_count=('litto','mean'),n_quadrats=('litto','size'))
assert snail_sites.n_quadrats.eq(3).all()
snail_sites.to_csv(OUT/'snail_site_times.csv',index=False)
pig = pd.read_excel(CACHE/'data_Pigments.xlsx')
endpoints = ['CarChl_ratio','Scytonemin']
assert len(pig)==24 and not pig.duplicated(['LightSite','rep']).any()
pig_sites = pig.groupby(['Light','LightSite','Herb'],as_index=False)[endpoints].mean()
assert pig.groupby(['Light','LightSite','Herb']).size().eq(3).all()
pig_sites.to_csv(OUT/'pigment_site_means.csv',index=False)
# Six variables per arm are NOT six independent lighting replicates.
# The four pigment variables are measured at the same two sites per arm.
cols = pd.MultiIndex.from_product([endpoints,['control','exclusion']])
arrays={}
for arm in ['dark','lit']:
    arrays[arm]=pig_sites[pig_sites.Light==arm].pivot(index='LightSite',columns='Herb',values=endpoints).reindex(columns=cols).to_numpy(float)
    assert arrays[arm].shape==(2,4) and np.isfinite(arrays[arm]).all()
mu={k:v.mean(axis=0) for k,v in arrays.items()}
theta=np.log(mu['lit']/mu['dark'])
V=sum(np.cov(x,rowvar=False,ddof=1)/2/np.outer(mu[k],mu[k]) for k,x in arrays.items())
ids=[f'Maggi2020_{e}_{g}' for e,g in cols]
def save(name,ids,est,cov,extra):
    se=np.sqrt(np.diag(cov));half=t.ppf(.975,1)*se
    rows=pd.DataFrame([dict(effect_id=id,lnRR=float(est[i]),sampling_variance=float(cov[i,i]),
                           ratio=float(np.exp(est[i])),approx_t1_lower=float(np.exp(est[i]-half[i])),
                           approx_t1_upper=float(np.exp(est[i]+half[i])),n_sites_lit=2,n_sites_dark=2,
                           main_catalogue=False,**extra[i]) for i,id in enumerate(ids)])
    rows.to_csv(OUT/f'{name}.csv',index=False)
    pd.DataFrame([dict(effect_id_i=ids[i],effect_id_j=ids[j],covariance=cov[i,j])
                  for i in range(len(ids)) for j in range(len(ids))]).to_csv(OUT/f'{name}_covariance.csv',index=False)
    assert np.linalg.eigvalsh(cov).min()>-1e-10
    return rows
save('pigment_light_effects',ids,theta,V,[dict(endpoint=e,grazer=g) for e,g in cols])
A=np.array([[-1.,1.,0.,0.],[0.,0.,-1.,1.]])
W=A@V@A.T
contrasts=save('pigment_exclusion_contrasts',[f'Maggi2020_{e}_exclusion_minus_control' for e in endpoints],A@theta,W,
              [dict(endpoint=e,estimand='lnRR_exclusion_minus_lnRR_grazed') for e in endpoints])
# Check covariance by independent linearized site scores, without pairing arms.
check=np.zeros((2,2))
for arm,x in arrays.items():
    influence=((x-mu[arm])/mu[arm])@A.T
    check+=np.cov(influence,rowvar=False,ddof=1)/2
assert np.allclose(W,check,rtol=1e-12,atol=1e-12)
last=snail_sites[snail_sites.time==snail_sites.time.max()]
assert last[last.herb=='e'].mean_count.eq(0).all()
lit=last[(last.light=='l2')&(last.herb=='c')].mean_count.to_numpy()
dark=last[(last.light=='l1')&(last.herb=='c')].mean_count.to_numpy()
assert len(lit)==len(dark)==2
v=lit.var(ddof=1)/2/lit.mean()**2+dark.var(ddof=1)/2/dark.mean()**2
save('final_snail_light_effect',['Maggi2020_day73_grazed_snails'],np.array([np.log(lit.mean()/dark.mean())]),np.array([[v]]),
     [dict(endpoint='snail_count_grazed_quadrats',day=73,mean_lit=lit.mean(),mean_dark=dark.mean())])
audit=dict(source_doi='10.5061/dryad.6m905qfw2',mirror='https://zenodo.org/records/5003041',
           files=manifest,snail_rows=216,snail_dates=9,pigment_rows=24,sites_per_arm=2,
           final_snail_day=73,excluded_snail_final_all_zero=True,
           pigment_effects=4,pigment_interaction_contrasts=2,snail_effects=1,
           covariance_influence_check=True,pigment_covariance_rank=int(np.linalg.matrix_rank(V)),
           cross_workbook_site_mapping='unresolved; no cross-endpoint covariance assumed',
           fluorescence_biomass_available=False,main_catalogue=False,
           unresolved=['Maggi2018 experimental-series overlap','Pigment-site to snail-site mapping',
                       'Fluorescence biomass and yield raw observations absent from recovered sheets'],
           interval_method='delta method with conservative t1; only two independent sites per arm',
           interpretation='Pigment ratio and scytonemin are not biomass, grazing or predation rates.')
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(contrasts[['effect_id','ratio','approx_t1_lower','approx_t1_upper']].to_string(index=False))
