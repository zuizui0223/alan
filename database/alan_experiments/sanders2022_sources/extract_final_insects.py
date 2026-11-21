"""Audit final field insect counts without inventing missing cages or covariance."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('source_dir', type=Path)
args = p.parse_args()
out = Path(__file__).resolve().parent
source = args.source_dir / 'Spectra_fieldexperiment.csv'
d = pd.read_csv(source)
columns = ['M.v', 'M.v.w', 'A.p', 'A.p.w', 'A.f', 'A.f.w', 'L.f', 'A.e', 'A.m']
assert len(d) == 850 and d.Cage.nunique() == 50
assert not d.duplicated(['Cage', 'Week']).any()
assert d.groupby('Cage').Week.nunique().eq(17).all()
assert d.groupby('Cage')[['Block', 'Treatment']].nunique().eq(1).all().all()
assert d['Abundance.factor'].isin([1, 2]).all()
assert d[columns].notna().all().all() and d[columns].ge(0).all().all()
f = d.loc[d.Week == d.Week.max()].copy()
assert len(f) == 50 and f.Date.nunique() == 1 and f.Date.iloc[0] == '4.10.17'
assert f['Abundance.factor'].eq(2).all()
for c in columns:
    f[c + '_expanded'] = f[c] * f['Abundance.factor']
f.to_csv(out / 'final_insect_cages.csv', index=False)
rows = []
for treatment, g in f.groupby('Treatment'):
    for c in columns:
        values = g[c + '_expanded']
        rows.append(dict(treatment=treatment, source_column=c, n_observed_cages=len(g),
                         blocks=';'.join(map(str, sorted(g.Block))),
                         expanded_mean=values.mean(), expanded_sd=values.std(ddof=1),
                         zero_cages=int(values.eq(0).sum()),
                         interpretation='column_specific_expanded_count_not_total_or_rate'))
summary = pd.DataFrame(rows)
summary.to_csv(out / 'final_insect_arm_summaries.csv', index=False)
contrasts = []
for treatment in sorted(set(f.Treatment) - {'C'}):
    lit = f.loc[f.Treatment == treatment].set_index('Block')
    dark = f.loc[f.Treatment == 'C'].set_index('Block')
    common = sorted(set(lit.index) & set(dark.index))
    for c in columns:
        lm = lit.loc[common, c + '_expanded'].mean()
        dm = dark.loc[common, c + '_expanded'].mean()
        finite = lm > 0 and dm > 0
        contrasts.append(dict(treatment=treatment, source_column=c,
            matched_blocks=';'.join(map(str, common)), n_matched_blocks=len(common),
            lit_mean_matched=lm, dark_mean_matched=dm,
            descriptive_ratio=lm/dm if dm > 0 else np.nan,
            descriptive_lnrr=np.log(lm/dm) if finite else np.nan,
            log_ratio_status='finite_point_only' if finite else 'zero_mean_not_log_transformable',
            pooling_status='hold_cohort_final_article_and_joint_covariance'))
contrasts = pd.DataFrame(contrasts)
contrasts.to_csv(out / 'final_insect_matched_points.csv', index=False)
audit = dict(source_file=source.name, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    source_rows=len(d), observed_cages=int(d.Cage.nunique()), final_week=17, final_date='2017-10-04',
    final_rows=len(f), expansion_factor=2, source_columns=columns,
    arm_summary_rows=len(summary), matched_contrast_rows=len(contrasts),
    finite_log_point_rows=int(contrasts.log_ratio_status.eq('finite_point_only').sum()),
    zero_mean_log_unavailable_rows=int(contrasts.log_ratio_status.ne('finite_point_only').sum()),
    observed_cages_by_treatment=f.groupby('Treatment').size().to_dict(),
    matched_n_by_treatment=contrasts.groupby('treatment').n_matched_blocks.first().to_dict(),
    missing_cages_not_filled=True, morph_columns_not_summed=True,
    new_pooled_effects=0, covariance_or_ci_estimated=False)
(out / 'final_insect_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
