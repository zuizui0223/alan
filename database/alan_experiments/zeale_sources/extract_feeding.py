"""Paired site-level descriptive contrasts; dark-first order remains confounded."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/zeale2018'
manifest = json.loads((CACHE / 'manifest.json').read_text())
for source in manifest:
    raw = (CACHE / source['file']).read_bytes()
    assert len(raw) == source['bytes']
    assert hashlib.sha256(raw).hexdigest() == source['sha256']
(OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
d = pd.read_csv(CACHE / 'bat_activity_data.csv')
keys = ['Site', 'Hedge side', 'Light type', 'Species']
assert len(d) == 320 and not d.duplicated(keys).any()
assert d.groupby('Site').size().eq(40).all() and d.Site.nunique() == 8
assert d['Number of bat passes'].sum() == 57558
g = d[(d.Species == 'Pipistrellus sp.') & (d['Hedge side'] == 'Lit')].copy()
assert len(g) == 40 and d['Number of feeding buzzes'].notna().sum() == 40
assert g['Number of feeding buzzes'].notna().all()
assert g['Number of bat passes'].gt(0).all()
g['exact_buzz_per_pass'] = g['Number of feeding buzzes'] / g['Number of bat passes']
assert (g.exact_buzz_per_pass-g['Feeding buzz ratio']).abs().le(0.00501).all()
g.to_csv(OUT / 'site_feeding_observations.csv', index=False)
rows, residuals, labels = [], [], []
endpoints = {'passes': 'Number of bat passes', 'buzzes': 'Number of feeding buzzes',
             'buzz_per_pass': 'exact_buzz_per_pass'}
for endpoint, col in endpoints.items():
    w = g.pivot(index='Site', columns='Light type', values=col).sort_index()
    assert w.shape == (8, 5) and not w.isna().any().any()
    dark = w.Dark.to_numpy()
    for light in ['Orange', 'White', 'Green', 'Red']:
        lit = w[light].to_numpy()
        yi = np.log(lit.mean()/dark.mean())
        influence = lit/lit.mean()-dark/dark.mean()
        vi = np.var(influence, ddof=1)/8
        half = t.ppf(.975, 7)*np.sqrt(vi)
        label = endpoint+'_'+light
        labels.append(label); residuals.append(influence)
        rows.append(dict(effect_id=label, endpoint=endpoint, light=light, paired_sites=8,
                         lit_mean=lit.mean(), dark_mean=dark.mean(), yi=yi, vi=vi,
                         ratio=np.exp(yi), lower=np.exp(yi-half), upper=np.exp(yi+half),
                         interval='paired_delta_t7', main_catalogue=False))
cov = np.cov(np.array(residuals), ddof=1)/8
assert np.allclose(np.diag(cov), [r['vi'] for r in rows])
assert np.linalg.eigvalsh(cov).min() > -1e-10
pd.DataFrame(rows).to_csv(OUT / 'feeding_contrasts.csv', index=False)
pd.DataFrame(cov,index=labels,columns=labels).to_csv(OUT / 'feeding_covariance.csv')
audit = dict(activity_rows=len(d),sites=8,passes_total=int(d['Number of bat passes'].sum()),
             feeding_rows=len(g),buzz_total=float(g['Number of feeding buzzes'].sum()),
             contrast_count=len(rows),covariance_rank=int(np.linalg.matrix_rank(cov)),
             max_original_ratio_rounding_error=float((g.exact_buzz_per_pass-g['Feeding buzz ratio']).abs().max()),
             dark_first_confounding_resolved=False,successful_captures_measured=False,
             adopted_main_effects=0)
(OUT / 'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(pd.DataFrame(rows)[['effect_id','ratio','lower','upper']].to_string(index=False))
