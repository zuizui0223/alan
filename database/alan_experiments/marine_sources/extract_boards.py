"""All measured taxa, two panels per board; exposure adjudication pending."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/marine2015'
metadata=json.loads((CACHE/'biostudies.json').read_text(encoding='utf-8'))
assert metadata['accno']=='S-EPMC4424621'
manifest=[]
for entry in metadata['section']['files']:
    raw=(CACHE/entry['path']).read_bytes();assert len(raw)==entry['size']
    manifest.append(dict(filename=entry['path'],bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
        url='https://www.ebi.ac.uk/biostudies/files/S-EPMC4424621/'+entry['path']))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
d=pd.read_csv(CACHE/'rsbl20150080supp3.csv')
assert len(d)==1692 and d['Unnamed: 0'].is_unique and d.notna().all().all()
keys=['Treatment','Board','Measure','Taxon']
assert d.groupby(keys).size().eq(2).all()
assert d.groupby('Board').Treatment.nunique().eq(1).all()
rows=[]
for k,g in d.groupby(keys,sort=True):
    value=g.Abundance.mean() if k[2]=='% cover' else g.Abundance.sum()
    assert k[2] in ['% cover','number of individuals']
    rows.append(dict(zip(keys,k),value=value))
b=pd.DataFrame(rows);b.to_csv(OUT/'board_measures.csv',index=False)
endpoints=sorted(set(zip(b.Measure,b.Taxon)))
assert len(endpoints)==47 and b.Board.nunique()==18
# No matching is inferred from board suffixes. Treatments were independent boards.
arms={arm:g.pivot(index='Board',columns=['Measure','Taxon'],values='value').reindex(columns=pd.MultiIndex.from_tuples(endpoints))
      for arm,g in b.groupby('Treatment')}
assert all(x.shape==(6,47) and x.notna().all().all() for x in arms.values())
order=['Control','Low','High'];means=np.concatenate([arms[a].mean().to_numpy() for a in order])
cov=np.zeros((141,141))
for i,a in enumerate(order):cov[i*47:(i+1)*47,i*47:(i+1)*47]=arms[a].cov().to_numpy()/6
effects=[];jac=[];ids=[];difference_jac=[]
for endpoint_index,(measure,taxon) in enumerate(endpoints):
    for arm,lux in [('Low',19),('High',30)]:
        li=order.index(arm)*47+endpoint_index;ci=endpoint_index
        a,c=means[li],means[ci];eid=f'Davies2015_E{endpoint_index+1:02d}_{lux}lux'
        jd=np.zeros(141);jd[li]=1;jd[ci]=-1;difference_jac.append(jd)
        result=dict(effect_id=eid,measure=measure,taxon=taxon,light_lux=lux,
            mean_lit=a,mean_dark=c,n_lit=6,n_dark=6,lnRR=np.nan,sampling_variance=np.nan,
            mean_difference=a-c,variance_difference=float(jd@cov@jd),
            lower=np.nan,upper=np.nan,df=np.nan,status='zero_arm_mean_no_lnRR')
        if a>0 and c>0:
            j=np.zeros(141);j[li]=1/a;j[ci]=-1/c
            variance=float(j@cov@j);estimate=np.log(a/c)
            va=cov[li,li]/a**2;vc=cov[ci,ci]/c**2
            df=variance**2/(va**2/5+vc**2/5) if variance>0 else np.nan
            half=t.ppf(.975,df)*np.sqrt(variance) if variance>0 else np.nan
            result.update(lnRR=estimate,sampling_variance=variance,lower=estimate-half,upper=estimate+half,
                df=df,status='finite_descriptive_exposure_hold' if variance>0 else 'zero_variance_hold')
            jac.append(j);ids.append(eid)
        effects.append(result)
e=pd.DataFrame(effects);e.to_csv(OUT/'all_taxon_contrasts.csv',index=False)
vd=np.asarray(difference_jac)@cov@np.asarray(difference_jac).T
assert np.allclose(np.diag(vd),e.variance_difference)
assert np.linalg.eigvalsh(vd).min()>-1e-8
pd.DataFrame(vd,index=e.effect_id,columns=e.effect_id).to_csv(OUT/'difference_covariance.csv',index_label='effect_id')
v=np.asarray(jac)@cov@np.asarray(jac).T
assert np.allclose(np.diag(v),e.set_index('effect_id').loc[ids].sampling_variance)
assert np.linalg.eigvalsh(v).min()>-1e-10
pd.DataFrame([dict(effect_id_1=a,effect_id_2=b,sampling_covariance=v[i,j])
 for i,a in enumerate(ids) for j,b in enumerate(ids)]).to_csv(OUT/'finite_covariance.csv',index=False)
audit=dict(raw_rows=len(d),boards=18,boards_per_arm=6,panels_per_board=2,
    measured_taxa=47,cover_taxa=10,count_taxa=37,all_contrasts=len(e),
    finite_lnRR=int(e.lnRR.notna().sum()),status_counts=e.status.value_counts().to_dict(),
    difference_contrasts=len(e),zero_variance_difference=int(e.variance_difference.eq(0).sum()),
    covariance_rank=int(np.linalg.matrix_rank(v)),catalogue_added=0,
    cautions=['All taxa retained; no significance or occurrence threshold used.',
              'All boards retained, including author-discussed outlier.',
              'No fabricated panel IDs; two panel measurements aggregate within each board.',
              'Methods dawn-on dusk-off wording conflicts with nocturnal description.',
              'Cover and numerical counts are separate scales; no global pooled effect.',
              'Zero treatment or control means are not replaced with pseudocounts.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
