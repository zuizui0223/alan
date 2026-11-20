"""Conditional cumulative oviposition effects from eight matched lighting blocks.

Missing daily observations are never filled. This reproduces source cumulative
columns, not a claim of seven complete daily counts or a rerun of the GLMM.
"""
from pathlib import Path
import hashlib
import json
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
RAW = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/daufel2023'
metadata = json.loads((RAW / 'metadata.json').read_text(encoding='utf8'))
f = metadata['files'][0]
source = RAW / f['name']
content = source.read_bytes()
assert hashlib.md5(content).hexdigest() == f['computed_md5']
arrays, rows, missing, units = [], [], [], []
for sheet, reported in [('Egg rafts', 261), ('Single eggs', 197)]:
    raw = pd.read_excel(source, sheet_name=sheet)
    d = raw[pd.to_numeric(raw.Pool, errors='coerce').isin(range(1, 17)) &
            pd.to_numeric(raw.Pair, errors='coerce').isin(range(1, 9))].copy()
    assert len(d) == 16 and d.Pool.is_unique
    d['Pair'] = pd.to_numeric(d.Pair).astype(int)
    assert set(d.Pair) == set(range(1, 9))
    assert d.groupby('Pair').Treatment.apply(set).eq({'Control', 'Lit'}).all()
    dates = [c for c in d if isinstance(c, datetime)]
    assert len(dates) == 7
    assert d[dates].notna().sum(axis=1).eq(6).all()
    for _, pair in d.groupby('Pair'):
        assert np.array_equal(pair[dates].isna().iloc[0], pair[dates].isna().iloc[1])
    observed_sum = d[dates].sum(axis=1, min_count=1)
    assert np.allclose(observed_sum, d.Cumulative)
    for _, r in d.iterrows():
        units.append(dict(endpoint=sheet, pair=r.Pair, pool=int(r.Pool), arm=r.Treatment,
                          cumulative=r.Cumulative, observed_days=int(r[dates].notna().sum())))
        for date in dates:
            if pd.isna(r[date]):
                missing.append(dict(endpoint=sheet, pair=r.Pair, pool=int(r.Pool),
                                    arm=r.Treatment, date=date.date().isoformat()))
    pivot = d.pivot(index='Pair', columns='Treatment', values='Cumulative').sort_index()
    a, b = pivot.Lit.to_numpy(float), pivot.Control.to_numpy(float)
    ma, mb = a.mean(), b.mean()
    assert ma > 0 and mb > 0
    # Pairwise influence on log ratio of arithmetic means preserves covariance.
    influence = a / ma - b / mb
    arrays.append(influence)
    variance = influence.var(ddof=1) / 8
    effect = np.log(ma / mb)
    half = t.ppf(.975, 7) * np.sqrt(variance)
    rows.append(dict(endpoint=sheet, n_pairs=8, mean_light=ma, mean_dark=mb,
                     lnRR=effect, variance=variance, ratio=np.exp(effect),
                     ci_lower=np.exp(effect-half), ci_upper=np.exp(effect+half),
                     source_total=int(a.sum()+b.sum()), published_total=reported,
                     total_matches=int(a.sum()+b.sum()) == reported,
                     pairs_light_lower=int((a < b).sum()),
                     catalogue_adopted=False))
cov = np.cov(np.array(arrays), ddof=1) / 8
assert np.linalg.eigvalsh(cov).min() > -1e-12
assert np.allclose(np.diag(cov), [r['variance'] for r in rows])
pd.DataFrame(rows).to_csv(OUT / 'conditional_effects.csv', index=False)
pd.DataFrame(units).to_csv(OUT / 'source_units.csv', index=False)
pd.DataFrame(missing).to_csv(OUT / 'missing_daily_cells.csv', index=False)
pd.DataFrame(cov, index=[r['endpoint'] for r in rows],
             columns=[r['endpoint'] for r in rows]).to_csv(OUT / 'effect_covariance.csv')
audit = dict(url=f['download_url'], filename=f['name'], md5=f['computed_md5'],
             sha256=hashlib.sha256(content).hexdigest(), bytes=len(content),
             missing_daily_cells=len(missing), missing_filled=False,
             inference='paired_delta_t7_conditional_source_cumulative',
             published_poisson_model_reproduced=False, final_adoption=False)
(OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf8')
print(pd.DataFrame(rows).to_string(index=False))
