"""Three-programme exploratory synthesis; see leaf_damage_meta_plan.md."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'leaf_damage_synthesis'
OUT.mkdir(exist_ok=True)
CACHE = ROOT.parents[1].parent/'gbif_annual/alan_experiment_sources/hongkong2025'
files = [ROOT/'cieraad_sources/descriptive_effects.csv', ROOT/'cieraad_sources/effect_covariance.csv',
         ROOT/'beech_sources/installation_effects.csv', ROOT/'beech_sources/installation_covariance.csv',
         CACHE/'Herb_Psyc.csv', CACHE/'Herb_Aporosa.csv']
parts, matrices, programmes = [], [], []
for name, file, covfile, endpoint, n in [
    ('Cieraad2023',files[0],files[1],None,2),
    ('Ceppi2025',files[2],files[3],'all_leaf_mean_percent',2)]:
    d = pd.read_csv(file)
    if endpoint: d = d[d.endpoint.eq(endpoint)].copy()
    d['programme'] = name
    parts.append(d[['effect_id','programme','lnRR','sampling_variance']])
    matrices.append(pd.read_csv(covfile,index_col=0).loc[d.effect_id,d.effect_id].to_numpy())
    programmes.append((name,n))
hk_rows, scores, site_rows = [], [], []
for file, species in zip(files[4:], ['Psychotria','Aporosa_genus_only']):
    d = pd.read_csv(file)
    d = d[d.Sampling.eq(3)].copy()
    assert len(d)==14 and not d.duplicated(['Site_pair','Setup']).any()
    assert d.groupby('Site_pair').Date.nunique().eq(1).all()
    g = d.groupby(['Site','Setup']).Percentage.mean().unstack().sort_index()
    assert len(g)==4 and g.notna().all().all() and g.gt(0).all().all()
    lit,dark = g.lit.to_numpy(),g.unlit.to_numpy()
    score = lit/lit.mean()-dark/dark.mean()
    scores.append(score)
    hk_rows.append(dict(effect_id='Guo2025_'+species,programme='Guo2025',
                        lnRR=np.log(lit.mean()/dark.mean()),sampling_variance=np.var(score,ddof=1)/4))
    for site,r in g.iterrows():site_rows.append(dict(species=species,site=site,mean_lit=r.lit,mean_dark=r.unlit))
parts.append(pd.DataFrame(hk_rows)); matrices.append(np.cov(np.array(scores),ddof=1)/4)
programmes.append(('Guo2025',4))
pd.DataFrame(site_rows).to_csv(OUT/'hongkong_site_means.csv',index=False)
contrasts=pd.concat(parts,ignore_index=True)
V=np.zeros((len(contrasts),len(contrasts))); A=np.zeros((3,len(contrasts)))
offset=0
for i,(part,cov) in enumerate(zip(parts,matrices)):
    m=len(part);V[offset:offset+m,offset:offset+m]=cov;A[i,offset:offset+m]=1/m;offset+=m
assert len(contrasts)==11 and np.allclose(np.diag(V),contrasts.sampling_variance)
assert np.linalg.eigvalsh(V).min()>-1e-10
y=A@contrasts.lnRR.to_numpy(); pv=A@V@A.T; v=np.diag(pv)
assert (v>0).all()
def nll(tau):
    vv=v+tau;w=1/vv;mu=np.sum(w*y)/w.sum()
    return .5*(np.log(vv).sum()+np.log(w.sum())+np.sum(w*(y-mu)**2))
bound=max(10.,10*np.var(y,ddof=1))
fit=minimize_scalar(nll,bounds=(0,bound),method='bounded',options={'xatol':1e-12})
assert fit.success and fit.x<bound*.99
tau=0. if nll(0)<=fit.fun else float(fit.x)
w=1/(v+tau);mu=float(np.sum(w*y)/w.sum());q=float(np.sum(w*(y-mu)**2)/2)
variance=float(max(1,q)/w.sum());radius=t.ppf(.975,2)*np.sqrt(variance)
pooled=dict(log_ratio=mu,ratio=np.exp(mu),lower=np.exp(mu-radius),upper=np.exp(mu+radius),
            tau2_REML=tau,q_over_df=q,mKH_variance=variance,df=2,programmes=3,
            dependent_contrasts=11,status='exploratory_not_final_adjudication')
pr=[]
for i,(name,n) in enumerate(programmes):
    rad=t.ppf(.975,n-1)*np.sqrt(v[i])
    pr.append(dict(programme=name,clusters=n,log_ratio=y[i],sampling_variance=v[i],
                   ratio=np.exp(y[i]),lower=np.exp(y[i]-rad),upper=np.exp(y[i]+rad),
                   meta_weight=w[i]/w.sum()))
contrasts['within_programme_weight']=contrasts.programme.map({name:1/len(part) for (name,_),part in zip(programmes,parts)})
contrasts.to_csv(OUT/'included_contrasts.csv',index=False)
pd.DataFrame(V,index=contrasts.effect_id,columns=contrasts.effect_id).to_csv(OUT/'sampling_covariance.csv')
pd.DataFrame(pv,index=[x[0] for x in programmes],columns=[x[0] for x in programmes]).to_csv(OUT/'programme_covariance.csv')
pd.DataFrame(pr).to_csv(OUT/'programme_estimates.csv',index=False)
pd.DataFrame([pooled]).to_csv(OUT/'pooled_estimate.csv',index=False)
audit=dict(input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
           covariance_ranks=[int(np.linalg.matrix_rank(c)) for c in matrices],
           final_eligibility=False,main_catalogue_changed=False,wetness_moderator_estimated=False,
           endpoint='terminal unconditional leaf-damage fractions in retained samples',
           estimand='equal log-contrast mean over each programme observed species-light menu',
           limitations=['2/2/4 clusters; uncertain delta sampling variances',
                        'Different species, spectra, dose menus, durations and leaf aggregation',
                        'Plant attrition and programme adjudication incomplete',
                        'No predation, chemical mediation or pollination mechanism estimated'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(pr).to_string(index=False));print(json.dumps(pooled,indent=2))
