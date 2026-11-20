# ALAN wet-habitat: is the strict matched light contrast geographically stable?

2026-11-20 | pre-outcome block influence and candidate/reference component check.
This is **not** a new causal test, a new independent botanical cohort, or
a 187-pair vegetation-adjusted analysis.

Existing original source SHA256 4d151698e0da5ee754326220d3bb817b141da00770299e7f8910366d8771e74f, commit 6a9a484.
Previously source-verified original matched samples: original 187/30 (failed
initial-water SMD gate), 0.5 robust-water caliper 133/26, 0.25 caliper 128/25.

Evaluate:
- Reproduce the identical original matched site IDs and numerical outcome
  decomposition. Report original and new pair memberships; a cohort effect
  difference may reflect **changed matching membership**, not mediation.
- Fit the same no-intercept covariate-adjusted light slope with baseline
  water, urbanization and climate and ANNUAL water-trend differences.
- Calculate candidate, reference, and candidate-minus-reference slope
  coefficients, with joint block-resampled 95% intervals (identical bootstrap
  draws preserve their covariance and algebraic decomposition).
- Leave exactly one 100-km source block out, refit the same coefficient,
  record range, sign persistence, and maximum influence. This is **not**
  holdout predictive validation of the exposure-to-flora relationship.
- Assess residual light exposure after linear projection on covariates,
  which diagnoses weak conditional identification but is not a p-value.

Scientific gates: reject changed source SHA, wrong 187/133/128 cohorts,
reused matched sites, insufficient rank, and failure of the candidate/
reference algebra. No invented coordinates or causal inference.
