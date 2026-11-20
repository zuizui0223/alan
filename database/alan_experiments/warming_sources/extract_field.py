"""Audit original allocation; extract only complete, unheated field stratum.

No treatment relabelling. Effects concern final aphid abundance, not predation.
Four neighbouring-row blocks are retained in joint delta-method covariance.
"""
from pathlib import Path
import hashlib
import json
import argparse
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--source', type=Path, default=REPO.parent / 'gbif_annual/alan_experiment_sources/warming_predation2017')
args = parser.parse_args()
source = args.source / 'Miller data ProcB.xlsx'
digest = hashlib.md5(source.read_bytes()).hexdigest()
assert digest == '85935aca65b83672b6b48fd3266d2fed'
dryad = json.loads((args.source / 'files.json').read_text())['_embedded']['stash:files'][0]
assert digest == dryad['digest'] and dryad['digestType'] == 'md5'
d = pd.read_excel(source, sheet_name='Field data')
assert len(d) == 48 and not d.duplicated(['Row', 'Cage']).any()
d['block'] = (d['Row'] + 1) // 2
d.to_csv(OUT / 'field_rows.csv', index=False)
expected = pd.MultiIndex.from_product([[1, 2, 3, 4], ['C7', 'CM', 'N'], ['N', 'Y'], ['N', 'Y']], names=['block', 'Pred.', 'Warm', 'light'])
counts = d.groupby(['block', 'Pred.', 'Warm', 'light']).size().reindex(expected, fill_value=0).rename('observed_n').reset_index()
counts.to_csv(OUT / 'allocation_audit.csv', index=False)
anomalies = counts[counts.observed_n != 1]
assert len(anomalies) == 2
# Restrict to the complete stratum, never fix the anomalous warmed labels.
u = d[d.Warm == 'N']
cols = pd.MultiIndex.from_product([['C7', 'CM', 'N'], ['N', 'Y']])
x = u.pivot(index='block', columns=['Pred.', 'light'], values='aphids').reindex(columns=cols).to_numpy(float)
assert x.shape == (4, 6) and np.isfinite(x).all()
mu = x.mean(axis=0)
cov_mean = np.cov(x, rowvar=False, ddof=1) / 4
J = np.zeros((3, 6))
effects = []
theta = np.log(mu[1::2] / mu[0::2])
for i, predator in enumerate(['C7', 'CM', 'N']):
    J[i, 2*i] = -1 / mu[2*i]
    J[i, 2*i+1] = 1 / mu[2*i+1]
    effects.append(dict(effect_id=f'Miller2017_unheated_{predator}', predator=predator,
                        endpoint='final_aphid_abundance', mean_dark=mu[2*i], mean_lit=mu[2*i+1],
                        sd_dark=x[:,2*i].std(ddof=1), sd_lit=x[:,2*i+1].std(ddof=1),
                        n_blocks=4, df=3, lnRR=theta[i], main_catalogue=False))
V = J @ cov_mean @ J.T
A = np.array([[1., 0., -1.], [0., 1., -1.]])
contrast = A @ theta
W = A @ V @ A.T
def save(rows, estimates, cov, stem):
    frame = pd.DataFrame(rows)
    frame['sampling_variance'] = np.diag(cov)
    half = t.ppf(.975, 3) * np.sqrt(np.diag(cov))
    frame['ratio'] = np.exp(estimates)
    frame['approx_t3_lower'] = np.exp(estimates - half)
    frame['approx_t3_upper'] = np.exp(estimates + half)
    frame.to_csv(OUT / f'{stem}.csv', index=False)
    ids = frame.effect_id.to_list()
    pd.DataFrame([dict(effect_id_i=ids[i], effect_id_j=ids[j], covariance=cov[i,j])
                  for i in range(len(ids)) for j in range(len(ids))]).to_csv(OUT / f'{stem}_covariance.csv', index=False)
    assert np.linalg.eigvalsh(cov).min() > -1e-10
    return frame
save(effects, theta, V, 'unheated_light_effects')
rows = [dict(effect_id=f'Miller2017_unheated_{p}_minus_none', predator=p,
             ln_ratio_of_response_ratios=contrast[i], n_blocks=4, df=3,
             endpoint='final_aphid_abundance_effect_modification', main_catalogue=False)
        for i,p in enumerate(['C7', 'CM'])]
result = save(rows, contrast, W, 'unheated_predator_contrasts')
# Independent scalar paired formulas check the matrix implementation.
for i in range(3):
    a,b = x[:,2*i],x[:,2*i+1]
    scalar = (a.var(ddof=1)/a.mean()**2 + b.var(ddof=1)/b.mean()**2
              - 2*np.cov(a,b,ddof=1)[0,1]/(a.mean()*b.mean()))/4
    assert np.isclose(scalar,V[i,i],rtol=1e-12,atol=1e-12)
influence = (x-mu) @ (A @ J).T
assert np.allclose(np.cov(influence,rowvar=False,ddof=1)/4,W,rtol=1e-12,atol=1e-12)
audit = dict(source_md5=digest, dryad_digest_match=True, source_doi='10.5061/dryad.h7k26',
             mirror='https://zenodo.org/records/4962533', paper_doi='10.1098/rspb.2017.1195',
             raw_rows=48, unique_cages=48, block_rule='ceil(Row/2); neighbouring rows in Methods',
             anomalous_allocation_cells=anomalies.to_dict('records'),
             extracted_stratum='Warm N only', unheated_rows=24, n_blocks=4,
             finite_light_effects=3, dependent_predator_contrasts=2,
             scalar_variance_check=True, contrast_influence_covariance_check=True,
             heated_stratum_held=True, main_catalogue=False,
             limitations=['Four blocks; delta-method t3 intervals are approximate.',
                          'Unheated stratum only; no warming interaction estimate.',
                          'Final abundance is not observed predation rate.',
                          'No correction of allocation labels; adjudication pending.'])
(OUT / 'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(result.to_string(index=False))
