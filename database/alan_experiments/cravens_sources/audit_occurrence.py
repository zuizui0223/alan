"""Reconstruct sample-level order occurrence; not prey abundance or pooled ALAN effects."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/cravens2018'
meta=json.loads((CACHE/'zenodo.json').read_text(encoding='utf-8'))
manifest=[]
for f in meta['files']:
 b=(CACHE/f['key']).read_bytes();assert 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
 manifest.append(dict(file=f['key'],source=f['links']['self'],bytes=len(b),md5=f['checksum'],sha256=hashlib.sha256(b).hexdigest()))
d=pd.read_excel(CACHE/'Cravens_et_al_2017_data.xlsx',sheet_name='2016_all_data')
s=d[['sampleid','batsp','lit','date']].drop_duplicates()
assert s.sampleid.is_unique
assert len(s)==188
assert not d['order'].isna().any()
expected={('EPFU','L'):7,('EPFU','U'):14,('LABO','L'):35,('LABO','U'):39,('MYGR','L'):7,('MYGR','U'):9,('MYLU','L'):9,('MYLU','U'):29,('NYHU','L'):6,('NYHU','U'):16,('PESU','L'):6,('PESU','U'):11}
assert s.groupby(['batsp','lit']).size().to_dict()==expected
# One detection per sample/order. Read depth is never interpreted as prey count.
a=d[['sampleid','batsp','lit','order']].drop_duplicates()
counts=a.groupby(['batsp','lit','order']).size().rename('samples_with_order').reset_index()
counts=counts.merge(s.groupby(['batsp','lit']).size().rename('analysed_samples').reset_index())
counts['order_detections_total']=counts.groupby(['batsp','lit']).samples_with_order.transform('sum')
counts['sample_prevalence']=counts.samples_with_order/counts.analysed_samples
counts['relative_order_occurrence']=counts.samples_with_order/counts.order_detections_total
counts.to_csv(OUT/'sample_order_occurrence.csv',index=False)
m=counts[counts.order=='Lepidoptera'].copy()
z=m.pivot(index='batsp',columns='lit',values='relative_order_occurrence')
z['lit_minus_unlit']=z.L-z.U
z.to_csv(OUT/'moth_occurrence_descriptive.csv')
audit=dict(raw_rows=len(d),samples=len(s),dates=s.date.nunique(),table3_all12cells_match=True,main_effects_added=0,definition='sample-order presence normalized over all order detections',figure1_reproduction=False,site_and_date_covariance='unresolved_no_meta_variance',paper_direction_check='MYLU decreases under written sample-order formula unlike reported direction; do not silently replace denominator',unique_denovo=int(d.denovo.nunique()))
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(z.to_string())
