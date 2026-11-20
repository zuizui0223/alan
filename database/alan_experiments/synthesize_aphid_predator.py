"""Two-programme exploratory synthesis of predator modification of ALAN counts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'aphid_predator_synthesis'
OUT.mkdir(exist_ok=True)
files = [ROOT/'bennie_trophic_sources/joint_predator_contrasts.csv',
         ROOT/'bennie_trophic_sources/joint_predator_contrast_covariance.csv',
         ROOT/'warming_sources/unheated_predator_contrasts.csv',
         ROOT/'warming_sources/unheated_predator_contrasts_covariance.csv']
p = pd.read_csv(files[0])
p = p[p.endpoint.eq('aphid_count')].copy()
p = p.rename(columns={'estimate':'log_ratio'})
p['endpoint'] = 'aphid_count_final_first_report_not_programme_terminal'
m = pd.read_csv(files[2]).rename(columns={'ln_ratio_of_response_ratios':'log_ratio'})
m['programme_id'] = 'Miller2017_unheated'
m['timepoint'] = '9_days_after_introduction'
m['endpoint'] = 'final_aphid_count'
m['light'] = 'LED_colour_unreported'
columns = ['effect_id','programme_id','endpoint','timepoint','light','log_ratio','sampling_variance']
d = pd.concat([p[columns], m[columns]], ignore_index=True)
assert len(p) == len(m) == 2 and d.effect_id.is_unique
v = np.zeros((4,4))
for indices, file in [(range(2),files[1]),(range(2,4),files[3])]:
    c = pd.read_csv(file).pivot(index='effect_id_i', columns='effect_id_j', values='covariance')
    ids = d.iloc[list(indices)].effect_id
    v[np.ix_(list(indices),list(indices))] = c.loc[ids,ids].to_numpy()
assert np.allclose(v,v.T) and np.allclose(np.diag(v), d.sampling_variance)
assert np.linalg.eigvalsh(v).min() > -1e-10
a = np.array([[.5,.5,0,0],[0,0,.5,.5]])
y = a @ d.log_ratio.to_numpy()
cv = a @ v @ a.T
assert abs(cv[0,1]) < 1e-12
vs = np.diag(cv)
tau2 = max(0., ((y[0]-y[1])**2-vs.sum())/2)
def restricted_nll(variance):
    vv = vs+variance
    w = 1/vv
    mu = np.sum(w*y)/w.sum()
    return .5*(np.log(vv).sum()+np.log(w.sum())+np.sum(w*(y-mu)**2))
fit = minimize_scalar(restricted_nll, bounds=(0,max(10.,4*tau2)), method='bounded', options={'xatol':1e-12})
numeric_tau2 = 0. if restricted_nll(0.) <= fit.fun else float(fit.x)
assert np.isclose(numeric_tau2,tau2,atol=1e-7)
w = 1/(vs+tau2)
mu = float(np.sum(w*y)/w.sum())
q = float(np.sum(w*(y-mu)**2))
variance = float(max(1.,q)/w.sum())
radius = float(t.ppf(.975,1)*np.sqrt(variance))
result = dict(log_ratio=mu,ratio=float(np.exp(mu)),lower=float(np.exp(mu-radius)),
    upper=float(np.exp(mu+radius)),tau2_REML=tau2,q=q,mKH_variance=variance,df=1,
    programmes=2,dependent_contrasts=4,status='exploratory_not_final_adjudication',
    endpoint='predator_modification_of_final_aphid_count_illumination_response')
rows = []
for i,(programme,df) in enumerate([('Penryn2012_2016',5),('Miller2017_unheated',3)]):
    rr = float(t.ppf(.975,df)*np.sqrt(vs[i]))
    rows.append(dict(programme_id=programme,log_ratio=y[i],sampling_variance=vs[i],
        ratio=np.exp(y[i]),approx_t_lower=np.exp(y[i]-rr),approx_t_upper=np.exp(y[i]+rr),
        programme_interval_df=df,meta_weight_fraction=w[i]/w.sum()))
pd.DataFrame([result]).to_csv(OUT/'pooled_estimate.csv',index=False)
pd.DataFrame(rows).to_csv(OUT/'programme_estimates.csv',index=False)
d['within_programme_weight'] = .5
d.to_csv(OUT/'included_contrasts.csv',index=False)
pd.DataFrame(v,index=d.effect_id,columns=d.effect_id).to_csv(OUT/'sampling_covariance.csv')
pd.DataFrame(cv,index=[r['programme_id'] for r in rows],columns=[r['programme_id'] for r in rows]).to_csv(OUT/'programme_covariance.csv')
audit = dict(programmes=2,contrasts=4,master_catalogue_changed=False,final_eligibility=False,
    input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
    analytic_tau2=tau2,numerical_tau2=numeric_tau2,
    covariance_min_eigenvalue=float(np.linalg.eigvalsh(v).min()),
    aggregation='equal log-contrast weights within observed experimental menus',
    interval='modified Knapp-Hartung max(1,q), t1',
    nonpooled_endpoints=['slug_period_mean','plant_final_biomass','aquatic_insect_emergence','pigment_ratio'],
    limitations=['Two programmes only, differing lights, predators, observation methods and durations',
        'Penryn 2013 is final first-report observation, not programme-terminal follow-up',
        'Sampling counts versus whole-cage census require comparable treatment-relative detection',
        'Miller non-warmed subset only; warmed treatment allocation unresolved',
        'Final sampled counts or census, not consumption rates',
        'No wetness moderator, pollination result, or contemporaneous plant mediation'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
print(json.dumps(result,indent=2))
