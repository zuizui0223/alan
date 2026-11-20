"""Render existing estimates only; no model fitting or new effect selection."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DB = ROOT / 'database/alan_experiments'
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
inputs = {}

def read(folder, filename):
    path = DB / folder / filename
    inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)

rows = []
def add(stage, label, row, pooled=False, lo='lower', hi='upper'):
    rows.append(dict(stage=stage, label=label, pooled=pooled,
                     ratio=float(row['ratio']), lower=row.get(lo, np.nan),
                     upper=row.get(hi, np.nan)))

c = read('capture_origin_synthesis', 'programme_estimates.csv')
c = c[c.variant.eq('family_counts_excluding_Ephydridae')]
assert len(c) == 2
assert c.programme_id.tolist() == ['Parkinson2016', 'SevenPonds2022']
for label, (_, r) in zip(['Parkinson 2020', 'Parkinson 2023'], c.iterrows()):
    add('A', label + ' (point only)', r)
p = read('capture_origin_synthesis', 'pooled_estimates.csv')
add('A', 'Pooled: 2 programmes', p[p.variant.eq('family_counts_excluding_Ephydridae')].iloc[0], True)
a = read('aphid_predator_synthesis', 'programme_estimates.csv')
assert len(a) == 2
assert a.programme_id.tolist() == ['Penryn2012_2016', 'Miller2017_unheated']
for label, (_, r) in zip(['Penryn', 'Miller 2017 (unheated)'], a.iterrows()):
    add('B', label, r, lo='approx_t_lower', hi='approx_t_upper')
add('B', 'Pooled: 2 programmes', read('aphid_predator_synthesis', 'pooled_estimate.csv').iloc[0], True)
l = read('leaf_damage_synthesis', 'programme_estimates.csv')
assert len(l) == 3
assert l.programme.tolist() == ['Cieraad2023', 'Ceppi2025', 'Guo2025']
for label, (_, r) in zip(['Netherlands', 'Switzerland', 'Hong Kong'], l.iterrows()):
    add('C', label, r)
add('C', 'Pooled: 3 programmes', read('leaf_damage_synthesis', 'pooled_estimate.csv').iloc[0], True)
d = pd.DataFrame(rows)
assert (d.ratio > 0).all()
bounds = d.dropna(subset=['lower', 'upper'])
assert ((bounds.lower > 0) & (bounds.lower <= bounds.ratio) & (bounds.ratio <= bounds.upper)).all()
d.to_csv(OUT / 'stage_meta_plot_data.csv', index=False)

titles = ['A  Aquatic versus terrestrial insect capture',
          'B  Predator modification of aphid response',
          'C  Terminal leaf damage']
xlabels = ['(Lit / dark, aquatic) / (lit / dark, terrestrial)',
           '(Lit / dark, predators present) / (lit / dark, predators absent)',
           'Lit / dark']
fig, axes = plt.subplots(3, 1, figsize=(11, 10))
fig.subplots_adjust(left=.27, right=.72, top=.9, bottom=.14, hspace=.72)
for ax, stage, title, xlabel in zip(axes, ['A', 'B', 'C'], titles, xlabels):
    sub = d[d.stage.eq(stage)].reset_index(drop=True)
    for i, r in sub.iterrows():
        y = len(sub) - 1 - i
        color = '#a3442d' if r.pooled else '#315e7a'
        if pd.notna(r.lower):
            ax.plot([r.lower, r.upper], [y, y], color=color, lw=1.7)
            value = f'{r.ratio:.3g} [{r.lower:.3g}, {r.upper:.3g}]'
        else:
            value = f'{r.ratio:.3g} [not displayed]'
        ax.scatter(r.ratio, y, marker='D' if r.pooled else 'o', s=45, color=color, zorder=3)
        ax.text(1.03, y, value, transform=ax.get_yaxis_transform(), va='center', fontsize=9)
    ax.set_yticks(list(reversed(range(len(sub)))))
    ax.set_yticklabels(sub.label, fontsize=9)
    ax.set_ylim(-.6, len(sub)-.4)
    ax.set_xscale('log')
    ax.set_xlim(.001, 2000)
    ax.set_xticks([.001,.01,.1,1,10,100,1000])
    ax.set_xticklabels(['0.001','0.01','0.1','1','10','100','1000'])
    ax.axvline(1, color='.5', ls='--', lw=.8)
    ax.axhline(.5, color='.85', lw=.7)
    ax.set_title(title, loc='left', fontsize=11, pad=12)
    ax.set_xlabel(xlabel + ' (log scale)', fontsize=9)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
fig.suptitle('Stage-specific experimental synthesis: provisional results', fontsize=14, y=.97)
fig.text(.03,.04,
         'Diamonds: REML pooled estimates with modified Knapp-Hartung 95% intervals.\n'
         'Programme intervals in B-C: existing small-sample t approximations. A: programme point estimates only;\n'
         'source tables provide variances but no comparable small-sample interval. No interval was newly constructed.\n'
         'Different endpoints are not pooled together. These comparisons do not test wet versus dry habitats.', fontsize=9)
for suffix in ['png','svg','pdf']:
    fig.savefig(OUT / f'stage_meta.{suffix}', dpi=200)
svg = OUT / 'stage_meta.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n', encoding='utf-8')
plt.close(fig)
(OUT / 'stage_meta_manifest.json').write_text(json.dumps(dict(
    input_sha256=inputs, plotted_rows=len(d), pooled_rows=int(d.pooled.sum()),
    refit=False, status='exploratory_not_final_adjudication',
    capture_programme_intervals='not displayed; no new interval constructed',
    common_axis_limits=[.001,2000]), indent=2)+'\n', encoding='utf-8')
print('Rendered 10 rows, including 3 separate pooled estimates; no refitting.')
