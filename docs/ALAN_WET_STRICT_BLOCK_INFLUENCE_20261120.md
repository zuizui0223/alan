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


## Executed original-source audit: 2026-11-20

[Original-source workflow #44674513433 (SUCCESS)](https://github.com/zuizui0223/alan/actions/runs/44674513433) passed all four unit tests and replayed the original historical source. The signed result archive is [artifact #14602355397](https://github.com/zuizui0223/alan/actions/runs/44674513433/artifacts/14602355397), ZIP SHA256 9600c44c33410f95a599eb1199662244acb2094a4a7190ba8d79e31859e73d74. All source coefficients, 100-km block omission runs and matching membership checks are [committed as JSON](../results/all_species_forest/wet_falsification/wet_strict_block_influence_components_20261120.json).

| Audit | 187 original pairs | 133 stricter pairs | 128 strictest pairs |
|---|---:|---:|---:|
| Spatial blocks | 30 | 26 | 25 |
| Same pair IDs retained from original 187 | 187 | 115 | 105 |
| New pair ID combinations | 0 | 18 | 23 |
| Full conditional candidate−reference coefficient, pp/10DN | −3.498 | −1.186 | −1.067 |
| Candidate component | −2.316 | −0.775 | −1.452 |
| Reference component | +1.181 | +0.411 | −0.386 |
| Candidate−reference joint-block bootstrap 95% CI | [−7.100,+0.545] | [−4.214,+1.880] | [−4.042,+1.981] |
| Block omissions still yielding a negative slope | **30/30** | **26/26** | **25/25** |
| Leave-one-block-out slope range | −4.411 to −1.491 | −2.079 to −0.369 | −1.898 to −0.345 |

Noteworthy: for the original 187 pairs, removing a 100-km block containing **just one matched pair** changes the conditional coefficient from −3.498 to −1.491. Thus the *magnitude* is sensitive to a single observation even though the *sign* remains negative for every omission. In the stringent 133/128 pairs, each individual block omission also yields a negative point slope, but this does **not** establish statistical significance: block-bootstrap CIs include zero throughout.

The candidate and reference components show no stable single explanation. In the 133 matched pairs, their point estimates were −0.775 and +0.411 pp/10DN; in the 128 sample, −1.452 and −0.386. All component 95% bootstrap CIs include zero. Do not equate a negative candidate−reference coefficient with observed abundance loss or sole candidate declines.

The stricter matched cohorts are different subpopulations, **not** just the original 187 unchanged comparisons: 133 retained 115 exact original pairings and added 18 alternative pairings; 128 retained 105 and added 23. Hence the attenuation in point estimates between columns does not measure hydrological mediation.

**Scientific disposition:** robust negative direction under single-block deletion, weak magnitude precision, one influential original block, and no demonstrated unique nightlight effect. This is an exploratory geographic stability diagnostic on postselected source wet candidates; original ALAN-independent experimental or held-out plant evidence remains necessary.
