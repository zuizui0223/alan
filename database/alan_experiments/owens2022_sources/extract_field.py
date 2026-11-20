"""Extract field mating counts only; do not treat pre-observations as trials."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    args = parser.parse_args()
    out = Path(__file__).resolve().parent
    raw = pd.read_excel(args.source, sheet_name='field')
    raw['source_excel_row'] = raw.index + 2
    assert set(raw.sequence) == {'pre', 'post'}
    data = raw.loc[raw.sequence.eq('post')].copy()
    # The supplied code defines success from time_to; independently compare
    # availability of stage_two before adopting this binary coding.
    assert data.time_to.notna().equals(data.stage_two.notna())
    data['reported_mate_success'] = data.time_to.notna().astype(int)
    assert set(data.treatment) == {'dark', 'light'}
    assert len(data) == 34 and data.date.nunique() == 14
    cols = ['source_excel_row', 'date', 'treatment', 'observer', 'time_start',
            'time_end', 'stage_one', 'stage_two', 'time_to', 'reported_mate_success']
    data[cols].to_csv(out / 'field_trials.csv', index=False)
    counts = data.groupby(['date', 'treatment'], as_index=False).agg(
        trials=('reported_mate_success', 'size'),
        successes=('reported_mate_success', 'sum'))
    counts.to_csv(out / 'field_date_arm_counts.csv', index=False)
    arms = data.groupby('treatment', as_index=False).agg(
        trials=('reported_mate_success', 'size'),
        successes=('reported_mate_success', 'sum'), dates=('date', 'nunique'))
    arms['proportion'] = arms.successes / arms.trials
    arms.to_csv(out / 'field_arm_counts.csv', index=False)
    p = arms.set_index('treatment').proportion
    audit = {
        'doi': '10.1098/rsos.220468',
        'source_url': 'https://doi.org/10.6084/m9.figshare.20407632',
        'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'source_rows': len(raw), 'excluded_pre_rows': int(raw.sequence.eq('pre').sum()),
        'post_trials': len(data), 'dates': int(data.date.nunique()),
        'dates_with_both_arms': int((counts.groupby('date').treatment.nunique() == 2).sum()),
        'descriptive_unadjusted_risk_ratio': float(p['light'] / p['dark']),
        'descriptive_unadjusted_risk_difference': float(p['light'] - p['dark']),
        'variance': None, 'meta_pool_eligible': False,
        'limitation': 'Date clustering and individual reuse unresolved; raw proportions are not adjusted model coefficients.',
    }
    (out / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
