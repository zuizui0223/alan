"""Pilot lnRR calculator: independent arms; uncertainty at the allocation unit.

No silent zero correction, inferred pairing, or conversion of individual SD
to plot SD. Other designs remain in the database for design-specific extraction.
"""
import csv
import math
from pathlib import Path
from scipy.stats import t, norm


def sd(row, arm):
    n = float(row[f'n_{arm}'])
    kind = row['dispersion_type']
    value = float(row[f'dispersion_{arm}'])
    if kind == 'SD':
        result = value
    elif kind == 'SE':
        result = value * math.sqrt(n)
    elif kind == 'CI_HALF_WIDTH':
        level = float(row['ci_level'])
        if not 0 < level < 1:
            raise ValueError('invalid CI level')
        distribution = row['ci_distribution']
        if distribution == 't':
            df = float(row[f'ci_df_{arm}'])
            if df <= 0:
                raise ValueError('invalid CI df')
            quantile = t.ppf((1 + level) / 2, df)
        elif distribution == 'normal':
            quantile = norm.ppf((1 + level) / 2)
        else:
            raise ValueError('CI distribution must be explicit')
        result = value / quantile * math.sqrt(n)
    else:
        raise ValueError('unsupported dispersion; do not infer')
    if not math.isfinite(result) or result < 0:
        raise ValueError('invalid dispersion')
    return result


def compute(row):
    if row['design'] != 'independent_arms' or row['dispersion_unit'] != row['independent_unit']:
        raise ValueError('design-specific extraction required')
    means = [float(row[f'mean_{a}']) for a in ('lit', 'dark')]
    ns = [float(row[f'n_{a}']) for a in ('lit', 'dark')]
    if any(not math.isfinite(m) or m <= 0 for m in means):
        raise ValueError('lnRR requires positive means; no pseudocount')
    if any(not math.isfinite(n) or n < 2 or n != int(n) for n in ns):
        raise ValueError('at least two independent units per arm required')
    deviations = [sd(row, a) for a in ('lit', 'dark')]
    variance = sum(s*s/(n*m*m) for s,n,m in zip(deviations,ns,means))
    return dict(lnRR=math.log(means[0]/means[1]), sampling_variance=variance,
                estimator='uncorrected_lnRR_first_order_delta', status='computed')


def main():
    root = Path(__file__).resolve().parent
    with (root/'effects.csv').open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    output = []
    for row in rows:
        try:
            result = compute(row)
        except (ValueError, KeyError) as exc:
            result = dict(lnRR='', sampling_variance='', estimator='', status=str(exc))
        output.append({**row, **result})
    with (root/'effects_computed.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print(f'{len(output)} effects processed; no pooled meta-analysis fitted')


if __name__ == '__main__':
    main()
