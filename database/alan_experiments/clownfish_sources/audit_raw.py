"""Audit original indoor breeding workbook; no field meta effects generated."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/clownfish2019'
meta=json.loads((CACHE/'zenodo.json').read_text(encoding='utf-8'))
f=meta['files'][0];source=CACHE/f['key'];b=source.read_bytes()
assert 'md5:'+hashlib.md5(b).hexdigest()==f['checksum']
book=pd.read_excel(source,sheet_name=None)
rows=[]
for name,df in book.items():
 pair='pair id' if name=='Hatch rate' else 'pairID'
 arm='treatment' if name=='Spawning interval' else 'treat'
 assert not df.duplicated([pair,'clutch','period']).any()
 for (a,t),g in df.groupby([arm,'period']):
  rows.append(dict(sheet=name,treatment=a.lower(),period=t,records=len(g),observed_pairs=g[pair].nunique(),pair_ids=';'.join(sorted(g[pair].unique()))))
fert=book['Fertilisation'];hatch=book['Hatch rate']
assert np.all(fert['egg.total']==fert.unfert+fert.fert)
assert np.allclose(fert['prop.fert'],fert.fert/fert['egg.total'])
assert np.all(hatch['num.hatched']<=hatch['num.healthy'])
assert np.allclose(hatch['prop.hatched'],hatch['num.hatched']/hatch['num.healthy'])
assert hatch.loc[(hatch.treat=='ALAN')&(hatch.period=='impact'),'num.hatched'].eq(0).all()
pd.DataFrame(rows).to_csv(OUT/'observed_support.csv',index=False)
# All observed clutches averaged within pair. Descriptive only, not final-time effects.
means=hatch.groupby(['pair id','treat','period'],as_index=False).agg(clutches=('clutch','size'),mean_hatch_proportion=('prop.hatched','mean'))
means.to_csv(OUT/'pair_period_descriptive.csv',index=False)
manifest=dict(source=f['links']['self'],bytes=len(b),md5=f['checksum'],sha256=hashlib.sha256(b).hexdigest(),metadata='https://zenodo.org/records/2542303')
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
audit=dict(sheet_rows={k:len(v) for k,v in book.items()},all_proportion_checks_pass=True,main_effects_added=0,lnRR_zero_handling='undefined_no_pseudocount',missing_pairs='not_imputed_as_zero',scope='indoor_context_only',pair_period_means='observed_clutches_only_not_intention_to_treat_or_final_time',final_eligibility=False)
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
