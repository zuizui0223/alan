import pathlib,re,json,hashlib
import pymupdf,numpy as np,pandas as pd
from scipy.stats import t
ROOT=pathlib.Path(__file__).resolve().parents[3]
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/toad2020'
OUT=pathlib.Path(__file__).resolve().parent
pdf=pymupdf.open(CACHE/'supplement.pdf')
text=' '.join(pdf[i].get_text() for i in [3,4])
num=r'(\d+(?:\.\d+)?)'
pat=r'\s+'.join([num]*10+[r'(on|off)',num])
rows=re.findall(pat,text)
assert len(rows)==37,len(rows)
cols=['day','month','mean_svl_mm','site','temperature_C','rainfall_mm','wind_kmh','lunar_percent','ambient_lux','ambient_category','lighting','mean_gut_g_per_toad_overnight']
d=pd.DataFrame(rows,columns=cols)
for c in cols:
 if c!='lighting':d[c]=pd.to_numeric(d[c])
assert sorted(d.groupby('site').size())==[6,6,6,6,6,7]
assert not d.duplicated(['month','day','site']).any()
d['source']='SupplementaryTableS2 pp4-5; text-table extraction visually checked'
d.to_csv(OUT/'trials_primary.csv',index=False)
g=d.groupby(['site','lighting']).agg(mean_gut_g=('mean_gut_g_per_toad_overnight','mean'),n_nights=('day','count'),ambient_lux=('ambient_lux','first')).reset_index()
g.to_csv(OUT/'enclosure_treatment_means.csv',index=False)
w=g.pivot(index='site',columns='lighting',values='mean_gut_g');x,y=w.on,w.off;n=len(w)
mx,my=x.mean(),y.mean();cov=np.cov(x,y,ddof=1)[0,1];inf=(x-mx)/mx-(y-my)/my;v=inf.var(ddof=1)/n;lnrr=np.log(mx/my);crit=t.ppf(.975,n-1)
out=dict(study_doi='10.1038/s41598-020-63503-9',endpoint='overnight prey gut mass g per toad; enclosure-treatment average over trial nights',window='Jan14-Apr3 2018; inferred year from main methods Dec2017-Apr2018',n_pairs=n,mean_lit=mx,sd_lit=x.std(ddof=1),mean_control=my,sd_control=y.std(ddof=1),within_enclosure_covariance=cov,lnRR=lnrr,variance_lnRR_paired=v,RR=mx/my,RR_t95_low=np.exp(lnrr-crit*np.sqrt(v)),RR_t95_high=np.exp(lnrr+crit*np.sqrt(v)),mean_difference_g=mx-my,variance_mean_difference=(x-y).var(ddof=1)/n,status='descriptive paired enclosure means; temporal treatment confounding remains')
pd.DataFrame([out]).to_csv(OUT/'effect_paired.csv',index=False)
(OUT/'source_manifest.json').write_text(json.dumps(dict(source='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7162902/supplementaryFiles',file='supplement.pdf',sha256=hashlib.sha256((CACHE/'supplement.pdf').read_bytes()).hexdigest(),doi='10.1038/s41598-020-63503-9',digitized=0,source_location='SupplementaryTableS2 pp4-5'),indent=2),encoding='utf-8')
print(g.to_string(index=False));print(json.dumps(out,indent=2))
