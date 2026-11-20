"""Final-week field responses with randomized-block and shared-control covariance."""
from pathlib import Path
import zipfile
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/sanders2018'
with zipfile.ZipFile(CACHE/'data.zip') as z:
    raw=z.read('data/ALANIntensityPopulationdynamics.csv')
    d=pd.read_csv(z.open('data/ALANIntensityPopulationdynamics.csv'))
assert len(d)==432 and not d.duplicated(['Week','Cage']).any()
assert d.groupby('Cage')[['Block','Treat.lux','Light.lux']].nunique().eq(1).all().all()
assert d.groupby('Week').size().eq(48).all()
last=d[d.Week==d.Week.max()].copy()
assert last.Week.eq(9).all() and last['ab.factor'].eq(1).all()
assert last.groupby(['Block','Treat.lux']).size().eq(1).all()
aphids=['M.v','M.v.w','A.p','A.p.w','A.f','A.f.w','S.a','S.a.w']
parasitoids=['L.f','A.e','A.m','P.dor.m','P.dor.a','P.dor.s']
assert last[aphids+parasitoids].notna().all().all()
last['aphid_abundance']=last[aphids].sum(axis=1)*last['ab.factor']
last['parasitoid_mummy_abundance']=last[parasitoids].sum(axis=1)*last['ab.factor']
metrics=['aphid_abundance','parasitoid_mummy_abundance']
last[['Week','Cage','Block','Treat.lux','Light.lux','ab.factor']+metrics].to_csv(OUT/'final_units.csv',index=False)
block=last.pivot(index='Block',columns='Treat.lux',values=metrics)
assert block.shape==(6,16) and block.notna().all().all()
means=block.mean();cov=block.cov().to_numpy()/6
rows=[];jac=[]
for metric in metrics:
    for dose in ['0.1','1','5','10','20','50','100']:
        a=(metric,dose);b=(metric,'c')
        j=np.zeros(len(means));j[means.index.get_loc(a)]=1/means[a];j[means.index.get_loc(b)]=-1/means[b]
        estimate=np.log(means[a]/means[b]);variance=float(j@cov@j)
        half=t.ppf(.975,5)*np.sqrt(variance)
        rows.append(dict(effect_id=f'Sanders2018_{metric}_{dose}lux',study_id='Sanders2018',
            endpoint=metric,light_lux=float(dose),mean_lit=means[a],mean_dark=means[b],
            n_lit=6,n_dark=6,paired_units=6,lnRR=estimate,sampling_variance=variance,
            lower=estimate-half,upper=estimate+half,experiment_id='Sanders2016_field48',
            control_group_id='Sanders2016_dark6',analysis_role='final_week_shared_control_doses',
            interpretation='Final-week standing counts; mummies are not per-host attack probability.'))
        jac.append(j)
e=pd.DataFrame(rows);v=np.asarray(jac)@cov@np.asarray(jac).T
assert np.allclose(np.diag(v),e.sampling_variance)
assert np.linalg.eigvalsh(v).min()>-1e-10
e.to_csv(OUT/'final_effects.csv',index=False)
pd.DataFrame([dict(effect_id_1=a,effect_id_2=b,sampling_covariance=v[i,j])
 for i,a in enumerate(e.effect_id) for j,b in enumerate(e.effect_id)]).to_csv(OUT/'final_covariance.csv',index=False)
manifest=[]
for filename,url in [('primary.pdf','https://kevingaston.com/wp-content/uploads/2019/05/Sanders-et-al-2018-Low-levels-of-artificial-light-at-night-change-food-web-dynamics.pdf'),
 ('data.zip','https://data-package.ceh.ac.uk/data/d30168d3-6cbb-4d75-b73c-276e6083a1fe.zip')]:
    b=(CACHE/filename).read_bytes();manifest.append(dict(filename=filename,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
audit=dict(source_rows=432,final_rows=48,blocks=6,final_week=9,effects=14,
    raw_csv_sha256=hashlib.sha256(raw).hexdigest(),covariance_rank=int(np.linalg.matrix_rank(v)),
    all_final_sampling_factors_one=True,pooled=False,
    caveats=['Final-week contrasts differ from published whole-season model.',
             'Two endpoints and seven doses share six randomized blocks and dark controls.',
             'Mummy abundance is not a parasitism rate or live parasitoid count.',
             'Missing plant-biomass labels and separate greenhouse/laboratory experiments not extracted here.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(e[['endpoint','light_lux','lnRR']].to_string(index=False))
