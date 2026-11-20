"""Audit complete growth observations; do not infer independent repeated cages."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/gurholt2026/original.csv'
assert hashlib.md5(SOURCE.read_bytes()).hexdigest() == 'e666e41643f7357c22fbea8345d83574'
d = pd.read_csv(SOURCE)
g = d.dropna(subset=['Initial Mass', 'Day 7 Mass']).copy()
assert not g.duplicated(['Cohort', 'Individual']).any()
for col in ['Hatch Date', 'Caterpillar Move Date', 'Day 7 date']:
    g[col + '_parsed'] = pd.to_datetime(g[col], format='mixed', errors='raise')
g['days_from_move'] = (g['Day 7 date_parsed'] - g['Caterpillar Move Date_parsed']).dt.days
g['days_from_hatch'] = (g['Day 7 date_parsed'] - g['Hatch Date_parsed']).dt.days
g['stage_normalized'] = g['Stage at move'].str.strip().str.lower().fillna('unknown')
g['treatment_conflict'] = g.Treatment.ne(g.Cage.str[:2])
g['mass_increment'] = g['Day 7 Mass'] - g['Initial Mass']
g['growth_divisor6'] = g.mass_increment / 6
g['growth_divisor7'] = g.mass_increment / 7
cols = ['Cohort', 'Individual', 'Cage', 'Treatment', 'Species', 'Hatch Date',
        'Caterpillar Move Date', 'Day 7 date', 'Stage at move', 'stage_normalized',
        'days_from_move', 'days_from_hatch', 'treatment_conflict', 'Initial Mass',
        'Day 7 Mass', 'mass_increment', 'growth_divisor6', 'growth_divisor7']
g[cols].to_csv(OUT / 'growth_interval_rows.csv', index=False)
cages = g.groupby(['Cohort', 'Cage', 'Treatment']).agg(
    complete_individuals=('Individual', 'size'),
    growth_divisor6=('growth_divisor6', 'mean'),
    growth_divisor7=('growth_divisor7', 'mean')).reset_index()
cages.to_csv(OUT / 'growth_cage_means.csv', index=False)
# Algebra check on the SAME sample and aggregation. This does not certify
# independence across cohorts or produce an adopted meta-analysis effect.
checks = []
for cohort, group in cages.groupby('Cohort'):
    estimates = []
    for col in ['growth_divisor6', 'growth_divisor7']:
        arms = group.groupby('Treatment')[col].agg(['mean', 'var', 'size'])
        y = np.log(arms.loc['AL', 'mean'] / arms.loc['ML', 'mean'])
        v = (arms['var'] / (arms['size'] * arms['mean'] ** 2)).sum()
        estimates.append((y, v))
    assert np.allclose(estimates[0], estimates[1], rtol=1e-12, atol=1e-14)
    checks.append({'cohort': int(cohort), 'max_absolute_difference':
                   float(np.max(np.abs(np.array(estimates[0]) - estimates[1])))})
audit = {
    'complete_growth_rows': len(g),
    'move_to_day7_days': {str(k): int(v) for k, v in g.days_from_move.value_counts().items()},
    'hatch_to_day7_days': {str(k): int(v) for k, v in g.days_from_hatch.value_counts().items()},
    'stage_counts': {k: int(v) for k, v in g.stage_normalized.value_counts().items()},
    'treatment_conflicts_in_growth_subset': int(g.treatment_conflict.sum()),
    'nonpositive_mass_increments': int(g.mass_increment.le(0).sum()),
    'cohort_cage_means_not_independent_installations': len(cages),
    'fixed_divisor_invariance_checks': checks,
    'initial_weighing_date_verified': False,
    'physical_cage_reuse_resolved': False,
    'adopted_effects': 0,
}
(OUT / 'growth_interval_audit.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
