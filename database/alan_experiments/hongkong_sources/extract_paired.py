"""Extract primary Figshare v1 plot-level endpoints; no pooled meta-analysis."""
import pathlib, json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import t
ROOT=pathlib.Path(__file__).resolve().parents[3]
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/hongkong2025'
OUT=pathlib.Path(__file__).resolve().parent
specs=[('spider_abundance','Nephila_Count.csv','Count'),('psychotria_herbivory','Herb_Psyc.csv','Percentage'),('aporosa_herbivory','Herb_Aporosa.csv','Percentage')]
rows=[];units=[];influences={}
for endpoint,filename,column in specs:
 d=pd.read_csv(CACHE/filename)
 assert not d.duplicated(['Plot','Date']).any()
 if endpoint=='spider_abundance':
  assert d.groupby('Setup')['Count'].sum().to_dict()=={'lit':201,'unlit':216}
  d=d[d.Date==d.groupby('Site_pair').Date.transform('max')]
 else:d=d[d.Sampling==3]
 assert len(d)==14 and d.Site_pair.nunique()==7
 assert d.groupby('Site_pair').Date.nunique().eq(1).all()
 w=d.pivot(index='Site_pair',columns='Setup',values=column).sort_index()
 x,y=w.lit,w.unlit;n=len(w);mx,my=x.mean(),y.mean()
 cov=np.cov(x,y,ddof=1)[0,1]
 inf=(x-mx)/mx-(y-my)/my
 v=inf.var(ddof=1)/n
 effect=np.log(mx/my);se=np.sqrt(v);crit=t.ppf(.975,n-1)
 rows.append(dict(endpoint=endpoint,source_file=filename,selection='last_same_date_per_pair' if endpoint=='spider_abundance' else 'round3',n_pairs=n,mean_lit=mx,sd_lit=x.std(ddof=1),mean_dark=my,sd_dark=y.std(ddof=1),within_pair_covariance=cov,lnRR=effect,variance_lnRR_paired=v,RR=np.exp(effect),RR_t95_low=np.exp(effect-crit*se),RR_t95_high=np.exp(effect+crit*se),mean_difference=(x-y).mean(),variance_mean_difference=(x-y).var(ddof=1)/n))
 influences[endpoint]=inf
 for _,r in d.iterrows():units.append(dict(endpoint=endpoint,pair=r.Site_pair,plot=r.Plot,date=r.Date,treatment=r.Setup,value=r[column]))
pd.DataFrame(rows).to_csv(OUT/'final_effects_paired.csv',index=False)
pd.DataFrame(units).to_csv(OUT/'final_plot_units.csv',index=False)
I=pd.DataFrame(influences)
(I.cov()/len(I)).to_csv(OUT/'final_effect_covariance.csv')
manifest=[]
for filename in [s[1] for s in specs]+['Spider_herbivory_to_ALAN_subtropical.zip','figshare_metadata.json']:
 f=CACHE/filename;manifest.append(dict(filename=filename,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size,source_doi='10.6084/m9.figshare.27273918.v1'))
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
