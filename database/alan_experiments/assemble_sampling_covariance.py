"""Assemble extracted sampling covariance with explicit programme dependence.

Not a meta-analysis fit. Singular aggregate/subset blocks are reported, not
silently regularized. Cross-programme sampling independence is an assumption.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'effect_catalogue.csv')
assert d.effect_id.is_unique
V=pd.DataFrame(0.0,index=d.effect_id,columns=d.effect_id)
audits=[]
for programme,g in d.groupby('programme_id',sort=False):
    assert g.covariance_table.nunique()==1,programme
    path=g.covariance_table.iloc[0]
    c=pd.read_csv(ROOT/path)
    if 'effect_id_1' in c:
        value='sampling_covariance' if 'sampling_covariance' in c else 'covariance'
        block=c.pivot(index='effect_id_1',columns='effect_id_2',values=value)
    elif 'effect_id_i' in c:
        block=c.pivot(index='effect_id_i',columns='effect_id_j',values='covariance')
    else:
        block=c.set_index(c.columns[0])
        lookup=dict(zip(g.endpoint,g.effect_id))
        assert set(block.index)==set(lookup) and set(block.columns)==set(lookup)
        block=block.rename(index=lookup,columns=lookup)
    ids=g.effect_id.tolist()
    block=block.loc[ids,ids]
    assert np.isfinite(block.to_numpy()).all()
    assert np.allclose(block,block.T,rtol=1e-8,atol=1e-12)
    assert np.allclose(np.diag(block),g.sampling_variance,rtol=1e-8,atol=1e-12)
    eigen=np.linalg.eigvalsh(block)
    assert eigen.min()>-1e-10,programme
    V.loc[ids,ids]=block.to_numpy()
    audits.append(dict(programme_id=programme,effects=len(ids),rank=int(np.linalg.matrix_rank(block)),
        min_eigenvalue=float(eigen.min()),source=path,
        full_rank=bool(np.linalg.matrix_rank(block)==len(ids))))
V.to_csv(ROOT/'sampling_covariance_all_extracted.csv',index_label='effect_id')
pd.DataFrame(audits).to_csv(ROOT/'sampling_covariance_block_audit.csv',index=False)
summary=dict(effect_rows=len(d),programmes=len(audits),
    diagonal_agreement=True,symmetric=True,positive_semidefinite=True,
    singular_programmes=[x['programme_id'] for x in audits if not x['full_rank']],
    fitted=False,
    warnings=['Cross-programme sampling covariance set to zero conditional on current programme identity audit.',
              'Shared totals and guild subsets can be linearly dependent; select a defined endpoint before fitting.',
              'Supplementary time means must not enter beside the same endpoint final-time effect as independent observations.',
              'Full rank does not establish eligibility, independence or valid moisture moderator support.'])
(ROOT/'sampling_covariance_audit.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
