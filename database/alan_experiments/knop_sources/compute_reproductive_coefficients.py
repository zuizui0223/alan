"""Intervals for printed model coefficients, not raw-arm lnRR estimates."""
from pathlib import Path
import csv
from statistics import NormalDist

out = Path(__file__).resolve().parent
with (out / 'reproductive_model_coefficients.csv').open(encoding='utf-8', newline='') as f:
    rows = list(csv.DictReader(f))
z = NormalDist().inv_cdf(0.975)
for row in rows:
    beta, se = float(row['estimate']), float(row['se'])
    row.update(variance=se**2, wald_lower=beta-z*se, wald_upper=beta+z*se,
               interval_method='normal_Wald_from_rounded_published_estimate_and_SE',
               cross_endpoint_covariance='unavailable')
with (out / 'reproductive_model_intervals.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(row['endpoint'], row['estimate'], row['wald_lower'], row['wald_upper'])
