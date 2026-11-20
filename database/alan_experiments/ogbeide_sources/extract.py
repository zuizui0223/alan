"""Extract supplied week-4 emergence snapshot at independent flume level.

Not the experiment's week-6 endpoint; keep outside final-time main catalogue.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
# Cache is alongside the repository, not under its tracked database.
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/ogbeide2026'
d = pd.read_csv(CACHE / 'selected/emergence_export.csv')
assert len(d) == 192
assert not d.duplicated(['week', 'flume', 'trap', 'day', 'Family']).any()
assert d.groupby('flume').treat.nunique().eq(1).all()
assert d.Family.eq('chiro').all()
assert np.allclose(d.abund_rate, d.abundance / (d.trap_day * d.area))
assert d.groupby(['week', 'flume']).size().eq(6).all()
g = d.groupby(['week', 'treat', 'flume'], observed=True).agg(
    capture_sum=('abundance', 'sum'), mean_rate=('abund_rate', 'mean'),
    trap_records=('abund_rate', 'size')).reset_index()
g.to_csv(OUT / 'flume_emergence_rates.csv', index=False)
last = g[g.week.eq('Week 4')]
rows = []
for name, lit, dark in [('no_crayfish', 'alan', 'control'),
                        ('with_crayfish', 'alan+cray', 'cray')]:
    a = last[last.treat.eq(lit)].mean_rate.to_numpy()
    b = last[last.treat.eq(dark)].mean_rate.to_numpy()
    assert len(a) == len(b) == 4
    y = np.log(a.mean() / b.mean())
    va = a.var(ddof=1) / (4 * a.mean()**2)
    vb = b.var(ddof=1) / (4 * b.mean()**2)
    variance = va + vb
    df = variance**2 / (va**2/3 + vb**2/3)
    q = t.ppf(.975, df) * np.sqrt(variance)
    rows.append(dict(effect_id='Ogbeide2026_week4_'+name, endpoint='chironomid_emergence_rate',
                     mean_lit=a.mean(), mean_dark=b.mean(), n_lit=4, n_dark=4,
                     lnRR=y, sampling_variance=variance, ratio=np.exp(y),
                     lower=np.exp(y-q), upper=np.exp(y+q), approximate_df=df,
                     analysis_role='week4_snapshot_not_final_week6'))
e = pd.DataFrame(rows)
e.to_csv(OUT / 'week4_light_effects.csv', index=False)
delta = e.lnRR.iloc[1] - e.lnRR.iloc[0]
variance = e.sampling_variance.sum()
summary = dict(doi='10.1111/1365-2435.70335', data_doi='10.5281/zenodo.14516568',
               source_rows=len(d), independent_flumes=d.flume.nunique(),
               measured_weeks=sorted(d.week.unique()), final_experiment_week=6,
               light_by_crayfish_lnRR_difference=float(delta),
               interaction_sampling_variance=float(variance),
               catalogue_included=False, pooled=False,
               warning='Source trap_day equals 4 for both d1 and d2. Do not infer calendar coverage beyond supplied sampling periods. Means use equal source effort; 6 trap observations are subsamples per flume. No inference about adult moth predation or pollination.')
(OUT/'extraction_audit.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
manifest=[]
for name in ['emergence.rds','emergence_export.csv']:
    f=CACHE/'selected'/name
    manifest.append(dict(file=name,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
                         source='Zenodo archive' if name.endswith('.rds') else 'R readRDS export of original archive object'))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
print(e.to_string(index=False))
print(json.dumps(summary,indent=2))
