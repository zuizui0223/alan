"""Original field slug damage: site-level extraction, no laboratory data."""
from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
S=Path('../gbif_annual/alan_experiment_sources/slug2025')
O=Path(__file__).resolve().parent
meta=json.loads((S/'dataset.json').read_text())
provenance=[]
for f in meta['files']:
 p=S/f['filename']
 if p.exists():
  h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==f['content_details']['sha256_hash']
  provenance.append(dict(file=p.name,sha256=h,url=f['content_details']['download_url']))
pd.DataFrame(provenance).to_csv(O/'source_manifest.csv',index=False)
# These five pair identities do not change across sampling periods in author code.
pairs={'MARTEN':('MUSTELID','LIT'),'BADGER':('MUSTELID','DARK'),'HORNET':('HYMENO','LIT'),'ANTS':('HYMENO','DARK'),'QUAIL':('WALKINGBIRD','LIT'),'UPUPA':('WALKINGBIRD','DARK'),'TOAD':('AMPHIBIAN','LIT'),'NEWT':('AMPHIBIAN','DARK'),'EAGLE':('RAPTOR','LIT'),'PEREGRINE':('RAPTOR','DARK')}
lit={'QUAIL','DEER','EAGLE','GRASS SNAKE','HORNET','MARTEN','TOAD','WOODPECKER'}
def labels(d):
 d['arm']=d['S_STU.SI'].map(lambda s:'LIT' if s in lit else 'DARK')
 d['stable_pair']=d['S_STU.SI'].map(lambda s:pairs.get(s,(None,None))[0])
 return d
f=labels(pd.read_csv(S/'C_240124_TRAITS.csv'))
f['date']=pd.to_datetime(f.S_DATE)
f['source_row']=np.arange(2,len(f)+2)
for part,cols in [('leaf',[f'M_LEA_DAM{i}' for i in range(10)]),('flower',[f'M_FLO_DAM{i}' for i in range(1,6)])]:
 a=f[cols].apply(pd.to_numeric,errors='raise')
 f[part+'_n_observed']=a.notna().sum(axis=1)
 f[part+'_proportion_damaged']=(a>0).sum(axis=1)/a.notna().sum(axis=1).replace(0,np.nan)
 f[part+'_mean_damage_all_observed']=a.mean(axis=1)
 # Positive-organ conditional severity: no undamaged organ replaced by missing data in unconditional endpoint.
 f[part+'_severity_damaged_only']=a.where(a>0).mean(axis=1)
cols=[c for c in f if c.startswith(('leaf_','flower_'))]
f[['source_row','P_PL_LABEL','P_PLANT_SP','S_STU.SI','date','arm','stable_pair']+cols].to_csv(O/'wild_plant_rows.csv',index=False)
means=[c for c in cols if not c.endswith('n_observed')]
wild=f.groupby(['P_PLANT_SP','S_STU.SI','arm','date'],dropna=False)[means].mean().reset_index()
wild['n_plants']=f.groupby(['P_PLANT_SP','S_STU.SI','arm','date'],dropna=False).size().values
wild['stable_pair']=wild['S_STU.SI'].map(lambda s:pairs.get(s,(None,None))[0])
wild.to_csv(O/'wild_site_date_means.csv',index=False)
# Original potted periods, author threshold 0.1 and clipping negative changes.
d=labels(pd.read_csv(S/'2023_E_HERBIVORY_FINAL_QUANTIFICATION-CLEAN.csv'))
d['date']=pd.to_datetime(d.S_DATE,dayfirst=True,format='mixed')
d['source_row']=np.arange(2,len(d)+2)
assert not d.duplicated(['P_PL_LABEL','date','S_TIME.PER']).any()
for when in ['BEFORE','AFTER']:
 a=d[[f'LEA_DAM{i}_{when}' for i in range(10)]]
 d['n_'+when]=a.notna().sum(axis=1)
 d['prop_'+when]=(a>.1).sum(axis=1)/a.notna().sum(axis=1).replace(0,np.nan)
 d['severity_'+when]=a.where(a>.1).mean(axis=1).mask(a.mean(axis=1)==0,0)
d['proportion_increase_raw']=d.prop_AFTER-d.prop_BEFORE
d['proportion_increase_author']=d.proportion_increase_raw.clip(lower=0)
d['severity_increase_author']=(d.severity_AFTER-d.severity_BEFORE).clip(lower=0)
d.to_csv(O/'potted_period_rows.csv',index=False)
metrics=['proportion_increase_author','severity_increase_author']
p=d.groupby(['S_STU.SI','arm','date','S_TIME.PER'])[metrics].mean().reset_index()
p['n_pot_periods']=d.groupby(['S_STU.SI','arm','date','S_TIME.PER']).size().values
p['stable_pair']=p['S_STU.SI'].map(lambda s:pairs.get(s,(None,None))[0])
p.to_csv(O/'potted_site_date_means.csv',index=False)
effects=[]; supports=[]; covariances=[]
def extract(data,tag,metrics):
 # Latest date with both arms within source-verified stable pair; not last individual plant.
 selected=[]
 for pair,g in data.dropna(subset=['stable_pair']).groupby('stable_pair'):
  good=[date for date,x in g.groupby('date') if set(x.arm)=={'LIT','DARK'}]
  if good:
   x=g[g.date==max(good)].copy();x['extraction']=tag;selected.append(x)
 if not selected:return
 a=pd.concat(selected);supports.append(a)
 influences={}
 for metric in metrics:
  z=a.pivot(index='stable_pair',columns='arm',values=metric).dropna()
  n=len(z)
  if n<2:continue
  x,y=z.LIT,z.DARK;mx,my=x.mean(),y.mean();cov=np.cov(x,y,ddof=1)[0,1]
  v=x.var(ddof=1)/(n*mx**2)+y.var(ddof=1)/(n*my**2)-2*cov/(n*mx*my) if mx>0 and my>0 else np.nan
  effects.append(dict(extraction=tag,endpoint=metric,n_pairs=n,mean_lit=mx,sd_lit=x.std(ddof=1),mean_dark=my,sd_dark=y.std(ddof=1),paired_covariance=cov,lnRR=np.log(mx/my) if mx>0 and my>0 else np.nan,variance=v,dates=';'.join(sorted(a.date.dt.strftime('%Y-%m-%d').unique())),status='holdout_source_identity_and_selected_support_review',degenerate_variance=bool(np.isfinite(v) and v<=1e-12)))
  if mx>0 and my>0:influences[metric]=x/mx-y/my
 for m,x in influences.items():
  for k,y in influences.items():
   if x.index.equals(y.index):covariances.append(dict(extraction=tag,endpoint1=m,endpoint2=k,n_pairs=len(x),covariance=x.cov(y)/len(x)))
for species,g in wild.groupby('P_PLANT_SP'):extract(g,'wild_'+species,means)
for period,g in p.groupby('S_TIME.PER'):extract(g,'potted_'+period,metrics)
pd.DataFrame(effects).to_csv(O/'paired_final_observed_effects.csv',index=False)
pd.DataFrame(covariances).to_csv(O/'same_support_effect_covariance.csv',index=False)
pd.concat(supports).to_csv(O/'selected_pair_site_support.csv',index=False)
audit=dict(doi='10.1016/j.isci.2025.112770',data_doi=meta['doi']['id'],wild_rows=len(f),wild_unique_labels=f.P_PL_LABEL.nunique(),wild_duplicate_labels=f[f.P_PL_LABEL.duplicated(False)][['P_PL_LABEL','S_STU.SI']].to_dict('records'),potted_rows=len(d),potted_unique_labels=d.P_PL_LABEL.nunique(),potted_unique_label_date=d[['P_PL_LABEL','date']].drop_duplicates().shape[0],potted_negative_proportion_changes=int((d.proportion_increase_raw<0).sum()),potted_missing_leaf_values=int(d.filter(regex=r'^LEA_DAM').isna().sum().sum()),laboratory_excluded=True,unobserved_sites_not_zero=True,authors_clipping_preserved=True)
(O/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
print(json.dumps(audit,indent=2));print(pd.DataFrame(effects)[['extraction','endpoint','n_pairs','lnRR','variance']].to_string(index=False))
