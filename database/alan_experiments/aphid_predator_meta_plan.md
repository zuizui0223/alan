# Exploratory predator modification of final aphid counts

2026-11-20. Written after inspecting the study-specific estimates and before
fitting the pooled model. This is not a preregistration or final review.

Endpoint: log illumination response ratio with predators minus log illumination
response ratio without predators. A ratio below one indicates a smaller final
aphid-count illumination response in the predator treatment; it is not an
observed predation rate or proof of mediation.

Include the two final aphid contrasts in Penryn (amber and white) and the two
non-warmed final aphid contrasts in Miller2017 (C7 and CM). Retain both directions
and use equal weights within each programme. The resulting estimand averages
the observed experimental menus: different colours in one study, different
predator species in the other. It is not standardized to a common predator,
light spectrum or duration. The programmes, not four contrasts, are the two
independent units of the between-study model. Propagate full within-programme
covariance; do not discard shared dark controls or common no-predator arms.

Use REML with its exact two-study solution, independently checked by numerical
minimization. Interval: modified Knapp-Hartung max(1,q), t1. Do not fit wetness
moderators, publication-bias tests or prediction intervals. Preserve the
programme-specific contrasts and explicitly report intervals spanning one.

Penryn aphid sampling and Miller cage census differ. Miller's complete unheated
stratum is selected because the warmed allocation cannot be reconciled; it is
not evidence for a warming effect. Do not use the warmer stratum's apparently
strong effects to select a result. Source labels, block reconstruction and final
source-version checks remain visible in the original extraction notes.

Penryn's 2013-09-12 count is the final observation in the first report,
Bennie2015, not the terminal observation of the entire programme. Subsequent
2014-2015 aphid records were recovered with the follow-up paper. Retention of
the first report's endpoint follows the review's first-report rule and must
not be represented as longest-follow-up evidence. Comparable relative detection
between lit/dark and predator/no-predator treatments is required to interpret
the sampling-count and census contrasts as abundance effect modification.

Slug period averages, aquatic emergence, pigment ratios and plant biomass are
different endpoints. In particular, Penryn aphids and later plant biomass do
not provide a contemporaneous mediation test. They are excluded from this
aphid synthesis, not discarded from the wider review.

Sources: [Bennie2015](https://doi.org/10.1098/rstb.2014.0131),
[Miller2017](https://doi.org/10.1098/rspb.2017.1195).
Method: [Rover, Knapp & Friede 2015](https://doi.org/10.1186/s12874-015-0091-1).
