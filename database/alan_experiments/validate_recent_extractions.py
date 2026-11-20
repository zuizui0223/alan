"""Independent computational cross-check, not human dual extraction.

Rebuild Hey pot means and Lu daytime block means from original workbooks.
Use scalar covariance identities independently of the extraction scripts.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE.parents[2] / 'gbif_annual/alan_experiment_sources'
checks = []

hey = pd.read_excel(RAW/'hey2020/supplement1.xlsx', sheet_name='allocation&defense')
he = pd.read_csv(HERE/'hey_sources/hey2020_biomass_effects.csv')
hc = pd.read_csv(HERE/'hey_sources/hey2020_biomass_covariance.csv')
values = {}
for r in he.itertuples():
    arms = []
    for arm in ['Lit', 'Control']:
        sub = hey[(hey.Light == arm) & (hey.Water == r.water_treatment) & (hey.Density == r.density)]
        means = {}
        for plot, rows in sub.groupby('Plot'):
            valid = rows.dropna(subset=['DryWeightS', 'DryWeightR'])
            means[plot] = sum(valid.DryWeightS + valid.DryWeightR) / len(valid)
        arms.append(pd.Series(means).sort_index())
    x, y = arms
    assert len(x) == len(y) == 5
    assert np.isclose(np.log(x.mean()/y.mean()), r.lnRR)
    v = x.var()/5/x.mean()**2 + y.var()/5/y.mean()**2
    assert np.isclose(v, r.variance)
    values[r.effect_id] = arms
for r in hc.itertuples():
    left, right = values[r.effect_id_i], values[r.effect_id_j]
    cv = sum(a.cov(b)/len(a)/a.mean()/b.mean() for a,b in zip(left,right))
    assert np.isclose(cv, r.covariance)
checks.append(dict(source='Hey2020', effects=len(he), covariance_entries=len(hc),
                   missing_root_individuals=int(hey.DryWeightR.isna().sum())))

lu = pd.read_excel(RAW/'lu2026/Date.xlsx', sheet_name='visitation')
day = lu[lu.sampling_time == 'Diurnal'].copy()
assert len(day)==200
assert not day.duplicated(['treatment','block','time','repeat']).any()
assert day.groupby(['block','treatment']).size().eq(20).all()
day['three_orders_DHC'] = day.D_abundance + day.H_abundance + day.C_abundance
le = pd.read_csv(HERE/'lu_sources/paired_block_descriptive_effects.csv')
le = le[(le.diel=='Diurnal') & (le.metric=='source_row_mean')]
lc = pd.read_csv(HERE/'lu_sources/daytime_effect_covariance.csv',index_col=0)
influences = {}
for r in le.itertuples():
    v = day.groupby(['block','treatment'])[r.endpoint].mean().unstack()
    x,y = v.Lit,v.Unlit
    z = x/x.mean()-y/y.mean()
    assert np.isclose(np.log(x.mean()/y.mean()),r.lnRR)
    assert np.isclose(z.var()/5,r.paired_delta_variance)
    influences[r.endpoint] = z
for a,x in influences.items():
    for b,y in influences.items():
        assert np.isclose(x.cov(y)/5,lc.loc[a,b])
checks.append(dict(source='Lu2026', effects=len(le), covariance_entries=lc.size,
                   independent_pairs=5, daytime_sessions=200))
report = dict(status='passed', checks=checks,
    scope='independent computational cross-check of original workbook arithmetic and dependence',
    not_verified=['human dual extraction','source discrepancy resolution','full eligibility adjudication'])
(HERE/'recent_extraction_validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
