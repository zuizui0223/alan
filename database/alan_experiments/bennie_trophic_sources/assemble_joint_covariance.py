"""Join Penryn responses by physical unit; never treat papers as experiments.

Different response times remain separate estimands. Joint sampling covariance
is not a mediation/path model and does not establish a trophic cascade.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
DB=OUT.parent
aph=pd.read_csv(DB/'bennie_sources/final_original_units.csv')
aph=aph[aph.endpoint=='aphid_count_3min'][['plot','trophic','light','value']]
biom=pd.read_csv(DB/'bennie_sources/anic2022_same_mesocosm_crosswalk.csv')
biom=biom.rename(columns={'trophic_2013':'trophic','light_2013':'light','total_biomass_g':'value'})[['plot','trophic','light','value']]
slug=pd.read_csv(OUT/'slug_original_unit_period_means.csv')
slug=slug.rename(columns={'Plot':'plot','T.treatment':'trophic','L.treatment':'light','period_mean':'value'})[['plot','trophic','light','value']]
slug.trophic=slug.trophic.map({'Grazer':'T2','Predator':'T3'})
slug.light=slug.light.replace({'Control':'Dark'})
units=[];rows=[];components={}
for endpoint,d,timepoint in [('aphid_count',aph,'2013-09-12'),('plant_biomass',biom,'2016_final_harvest'),('slug_count',slug,'15_month_mean_2013_2014')]:
    d=d.copy();d['endpoint']=endpoint;d['timepoint']=timepoint;units.append(d)
    for trophic,g in d.groupby('trophic'):
        if endpoint=='aphid_count' and trophic=='T1':continue
        dark=g[g.light=='Dark'].set_index('plot').value
        assert len(dark)==6 and dark.notna().all()
        for light in ['Amber','White']:
            lit=g[g.light==light].set_index('plot').value
            assert len(lit)==6 and lit.notna().all() and min(lit.mean(),dark.mean())>0
            if endpoint=='aphid_count':id=f'Bennie2015_aphid_count_3min_{trophic}_{light}'
            elif endpoint=='plant_biomass':id=f"Anic2022_community_{trophic}_{ {'Amber':'L2','White':'L3'}[light]}"
            else:id=f"Bennie2018_slug_period_{ {'T2':'Grazer','T3':'Predator'}[trophic]}_{light}"
            components[id]={(trophic,light):lit/lit.mean(),(trophic,'Dark'):-dark/dark.mean()}
            rows.append(dict(effect_id=id,endpoint=endpoint,trophic=trophic,light=light,timepoint=timepoint,
                             lnRR=np.log(lit.mean()/dark.mean()),n_per_cell=6,programme_id='Penryn2012_2016',
                             role='secondary_period_mean' if endpoint=='slug_count' else 'existing_final_endpoint'))
edf=pd.DataFrame(rows);ids=edf.effect_id.tolist()
def cov(a,b):
    value=0.
    for cell in a.keys() & b.keys():
        x,y=a[cell].align(b[cell],join='inner')
        assert len(x)==len(y)==6
        value+=np.cov(x,y,ddof=1)[0,1]/6
    return value
V=np.array([[cov(components[a],components[b]) for b in ids] for a in ids])
edf['sampling_variance']=np.diag(V)
edf.to_csv(OUT/'joint_effects.csv',index=False)
pd.concat(units).to_csv(OUT/'joint_original_units.csv',index=False)
def save_cov(name,ids,V):
    pd.DataFrame([dict(effect_id_i=a,effect_id_j=b,covariance=V[i,j]) for i,a in enumerate(ids) for j,b in enumerate(ids)]).to_csv(OUT/name,index=False)
save_cov('joint_sampling_covariance.csv',ids,V)
assert np.allclose(V,V.T) and np.linalg.eigvalsh(V).min()>-1e-10
# Cross-check both already published covariance submatrices exactly.
old=pd.read_csv(DB/'bennie_sources/programme_covariance_including_anic2022.csv')
for r in old.itertuples():
    assert np.isclose(V[ids.index(r.effect_id_1),ids.index(r.effect_id_2)],r.sampling_covariance,atol=1e-12,rtol=1e-12)
old_slug=pd.read_csv(OUT/'slug_period_effects_covariance.csv')
for r in old_slug.itertuples():
    assert np.isclose(V[ids.index(r.effect_id_i),ids.index(r.effect_id_j)],r.covariance,atol=1e-12,rtol=1e-12)
contrasts=[];A=[]
for endpoint in ['aphid_count','slug_count','plant_biomass']:
    for light in ['Amber','White']:
        v=np.zeros(len(ids))
        for trophic,sign in [('T3',1),('T2',-1)]:
            row=edf[(edf.endpoint==endpoint)&(edf.light==light)&(edf.trophic==trophic)].iloc[0]
            v[ids.index(row.effect_id)]=sign
        A.append(v)
        contrasts.append(dict(effect_id=f'Penryn_{endpoint}_predator_{light}',endpoint=endpoint,light=light,
                              timepoint=row.timepoint,estimand='lnRR_predator_minus_lnRR_herbivore_only'))
A=np.array(A);theta=A@edf.lnRR.to_numpy();W=A@V@A.T
out=pd.DataFrame(contrasts);out['estimate']=theta;out['sampling_variance']=np.diag(W)
half=t.ppf(.975,5)*np.sqrt(np.diag(W))
out['ratio_of_response_ratios']=np.exp(theta)
out['approx_t5_lower']=np.exp(theta-half);out['approx_t5_upper']=np.exp(theta+half)
out['programme_id']='Penryn2012_2016';out['pooled']=False
out.to_csv(OUT/'joint_predator_contrasts.csv',index=False)
save_cov('joint_predator_contrast_covariance.csv',out.effect_id.tolist(),W)
# Independent centered-score Gram matrix on54 physical plots.
plots=sorted(biom['plot']);scores=np.zeros((len(plots),len(ids)))
for j,id in enumerate(ids):
    for values in components[id].values():
        for plot,val in (values-values.mean()).items():scores[plots.index(plot),j]=val/np.sqrt(6*5)
assert np.allclose(scores.T@scores,V,atol=1e-12,rtol=1e-12)
assert np.allclose((scores@A.T).T@(scores@A.T),W,atol=1e-12,rtol=1e-12)
audit=dict(physical_mesocosms=54,programme_count=1,study_labels=3,effects=14,
           joint_covariance_cells=196,predator_contrasts=6,contrast_covariance_cells=36,
           existing_100_cell_submatrix_verified=True,slug_16_cell_submatrix_verified=True,
           independent_centered_score_check=True,min_eigenvalue=float(np.linalg.eigvalsh(V).min()),
           rank=int(np.linalg.matrix_rank(V)),main_catalogue_unchanged=True,
           new_independent_experiments=0,pooled_estimates=0,
           limitations=['Different response periods; no concurrent chain or mediation estimate',
                        'Slug period mean remains supplementary, not final endpoint',
                        'Empirical covariance from six plots per cell; approximate t5 intervals',
                        'Intermediate aphid version discrepancies retained separately'])
(OUT/'joint_covariance_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(out[['endpoint','light','ratio_of_response_ratios','approx_t5_lower','approx_t5_upper']].to_string(index=False))
