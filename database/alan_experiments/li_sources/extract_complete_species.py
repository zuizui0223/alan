"""Audit publicly deposited preprint data; complete-design subset only.

Raw table is not redistributed (author README: data for peer review only).
All estimates remain descriptive and outside the main catalogue.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/li2026'
raw=(CACHE/'seedset.csv').read_bytes()
assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()=='b616aca5beeba1457549841477d61ca4f796a147'
d=pd.read_csv(CACHE/'seedset.csv')
counts=d.groupby(['species','Treatment','replicate']).size().rename('rows').reset_index()
counts.to_csv(OUT/'source_cell_counts.csv',index=False)
complete=[]
for species,g in d.groupby('species'):
    expected_rounds=3 if species=='ML' else 1
    cells=g.groupby(['round','Treatment','replicate']).size()
    if len(cells)==9*expected_rounds and cells.eq(2).all() and g[['Seeds','Fruit','pollinated']].notna().all().all():
        complete.append(species)
assert complete==['SN']
g=d[d.species.eq('SN')].copy()
assert g.id.is_unique and len(g)==18
assert g.pollinated.eq(g.Seeds.gt(0).astype(int)).all()
block=g.groupby(['replicate','Treatment'])[['Seeds','Fruit','pollinated']].mean().unstack('Treatment')
block.to_csv(OUT/'SN_block_means.csv')
means=block.mean(); cov=block.cov().to_numpy()/len(block)
effects=[];jac=[]
for metric in ['Seeds','Fruit','pollinated']:
    for treatment in ['L','P']:
        j=np.zeros(len(means));a=means[(metric,treatment)];b=means[(metric,'C')]
        if metric=='pollinated':
            value=a-b;j[means.index.get_loc((metric,treatment))]=1;j[means.index.get_loc((metric,'C'))]=-1;scale='risk_difference'
        else:
            assert a>0 and b>0
            value=np.log(a/b);j[means.index.get_loc((metric,treatment))]=1/a;j[means.index.get_loc((metric,'C'))]=-1/b;scale='lnRR'
        variance=float(j@cov@j)
        q=t.ppf(.975,2)*np.sqrt(max(variance,0))
        effects.append(dict(effect_id=f'Li2026_SN_{treatment}_{metric}',endpoint=metric,
            treatment=treatment,n_blocks=3,mean_lit=a,mean_dark=b,scale=scale,
            estimate=value,sampling_variance=variance,lower=value-q,upper=value+q,
            status='descriptive_preprint_complete_species_subset_not_catalogue'))
        jac.append(j)
e=pd.DataFrame(effects);e.to_csv(OUT/'SN_descriptive_effects.csv',index=False)
v=np.asarray(jac)@cov@np.asarray(jac).T
assert np.allclose(np.diag(v),e.sampling_variance) and np.linalg.eigvalsh(v).min()>-1e-10
pd.DataFrame(v,index=e.effect_id,columns=e.effect_id).to_csv(OUT/'SN_effect_covariance.csv',index_label='effect_id')
audit=dict(preprint_doi='10.64898/2026.07.10.737511',source_blob='b616aca5beeba1457549841477d61ca4f796a147',
    source_sha256=hashlib.sha256(raw).hexdigest(),source_rows=len(d),paper_planned_plants=90,
    rows_by_species=d.species.value_counts().to_dict(),missing_seed_outcomes=int(d.Seeds.isna().sum()),
    duplicate_label_rows=int(d[d.Label.notna()].duplicated(['species','round','Label'],keep=False).sum()),
    complete_species_codes=complete,effect_rows=len(e),independent_blocks=3,
    covariance_rank=int(np.linalg.matrix_rank(v)),catalogue_included=False,
    holds=['Preprint; species codes matched by paper species initials.',
           'ML and ZC row/label discrepancies unresolved; no relabelling or zero imputation.',
           'SN complete cells do not resolve entire study audit; three blocks and approximate intervals.',
           'Lighting arrangement also changes lamp count; no direct moth or day/night pollination measurement.',
           'Published plant-level model uncertainty not substituted for block-level uncertainty.',
           'Raw data not redistributed; public author README labels data for peer review only.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(e.to_string(index=False))
