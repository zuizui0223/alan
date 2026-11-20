"""Conditional mass gain by cohort; reused cages are never pooled as studies."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/gurholt2026/original.csv'
assert hashlib.md5(SOURCE.read_bytes()).hexdigest() == 'e666e41643f7357c22fbea8345d83574'
d = pd.read_csv(SOURCE)
assert not d.duplicated(['Cohort', 'Individual']).any()
g = d.dropna(subset=['Initial Mass', 'Day 7 Mass']).copy()
assert g.Treatment.eq(g.Cage.str[:2]).all()
g['mass_gain'] = g['Day 7 Mass'] - g['Initial Mass']
assert (g.mass_gain > 0).all()
g['stage'] = g['Stage at move'].str.strip().str.lower().fillna('unknown')
cages = g.groupby(['Cohort', 'Cage', 'Treatment']).agg(
    complete_individuals=('Individual', 'size'), mean_gain=('mass_gain', 'mean'),
    host_inc=('Species', lambda s: int(s.eq('inc').sum())),
    host_syr=('Species', lambda s: int(s.eq('syr').sum())),
    first_instar=('stage', lambda s: int(s.eq('first instar').sum())),
    second_instar=('stage', lambda s: int(s.eq('second instar').sum())),
    unknown_stage=('stage', lambda s: int(s.eq('unknown').sum()))).reset_index()
cages.to_csv(OUT/'mass_gain_cage_units.csv', index=False)
out = []
for cohort, block in cages.groupby('Cohort'):
    arms = block.groupby('Treatment').mean_gain.agg(['mean', 'var', 'size'])
    expected = 3 if cohort == 1 else 4
    assert set(arms.index) == {'AL', 'ML'} and arms['size'].eq(expected).all()
    means = arms['mean']
    parts = arms['var']/(arms['size']*means**2)
    y = float(np.log(means.AL/means.ML))
    v = float(parts.sum())
    df = float(v**2 / (parts**2/(arms['size']-1)).sum())
    half = t.ppf(.975, df)*np.sqrt(v)
    out.append(dict(cohort=int(cohort), endpoint='recorded_mass_gain_conditional_on_two_weighings',
        mean_lit=float(means.AL), mean_dark=float(means.ML), n_cages_per_arm=expected,
        complete_individuals=int(block.complete_individuals.sum()),
        lnRR=y, sampling_variance=v, welch_df=df, ratio=float(np.exp(y)),
        approx95_lower=float(np.exp(y-half)), approx95_upper=float(np.exp(y+half)),
        role='descriptive_within_cohort_not_independent_programmes', pooled=False))
pd.DataFrame(out).to_csv(OUT/'cohort_mass_gain_effects.csv', index=False)
audit = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), complete_rows=len(g),
    cohort_cage_rows=len(cages), supported_cohorts=sorted(map(int,g.Cohort.unique())),
    unsupported_cohorts=sorted(set(map(int,d.Cohort.unique()))-set(map(int,g.Cohort.unique()))),
    reused_cage_labels=len(set(cages.Cage)), cross_cohort_covariance='not_assumed_zero_not_estimated',
    grouping='equal_weight_cage_means_not_individuals', day_divisor='none',
    caveats=['Conditional on observed initial and followup masses; missing is not death or zero',
             'Different stages and host composition retained; no covariate-adjusted causal claim',
             'Alternating light allocation; repeated physical cages and date errors unresolved',
             'Approximate delta/Welch intervals describe cage variation within each cohort',
             'No survival, adult performance, pollination or population growth inference'],
    main_catalogue_rows_added=0, independent_programmes_added=0)
(OUT/'cohort_mass_gain_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(out)[['cohort','complete_individuals','ratio','approx95_lower','approx95_upper']].to_string(index=False))
