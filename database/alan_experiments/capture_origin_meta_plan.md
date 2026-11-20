# Exploratory capture-origin synthesis specification

Written before fitting the pooled model on 2026-11-20, after inspecting the
individual estimates. This is an exploratory analysis, not a preregistration.

Question: among the two currently extracted freshwater-shoreline experimental
programmes, how large is the relative illumination response of aquatic-origin
versus terrestrial-origin captures?

Use every finite contrast in `capture_origin_contrasts.csv`; do not select by
direction, significance, colour, or placement. Within each programme take the
equal-weight average of the log response-ratio differences. In SevenPonds2022
this averages four colours and two shoreline placements. This deliberately
defines a treatment-menu average; it does not estimate a common white-light
effect or a wet-versus-dry habitat contrast. Transform the full supplied sampling
covariance by the same weights. The two programmes, not nine contrasts, are the
units of the between-programme model. The different menus and single-night
capture windows limit generalization.

Fit an intercept-only random-effects model by REML. With exactly two independent
programme estimates, use the closed-form REML variance estimate and verify it
against numerical restricted-likelihood minimization. Report a modified
Knapp-Hartung interval using max(1, q) and t with k-1 = 1 degrees of freedom;
do not substitute a narrower normal interval when the modified interval crosses
one. Do not fit wetness moderators, publication-bias tests, or a prediction
interval with these two programmes. A boundary estimate of heterogeneity does
not establish homogeneity.

Source/design adjudication remains provisional: Parkinson2016 fish-allocation
identities are unavailable; SevenPonds2022 has only three stations per cell,
possible between-station attraction, and coarse taxon-origin assignments.
Do not mark the master catalogue as finally accepted. Save this as a distinct
exploratory supplementary synthesis and preserve the programme-specific results.

The pooled endpoint is a ratio of capture response ratios. It is neither
predation nor pollination, and does not demonstrate a mechanism of plant decline.

Method reference: [Rover, Knapp & Friede 2015](https://doi.org/10.1186/s12874-015-0091-1).
Primary sources: [Parkinson2020](https://doi.org/10.1371/journal.pone.0240138)
and [Parkinson2023](https://doi.org/10.1098/rstb.2022.0361).

## Amendment after source review, before successful pooled fit

The independent source audit found Ephydridae assigned to aquatic origin in
2020 and terrestrial origin in 2023, plus four 2020 mesocosms whose aquatic
summary totals differ from the sum of family columns. Retain original outputs;
create a family-count version excluding the conflicting family from both studies.
This is one documented conflict among exact-name matches, not full taxonomic
harmonization. No unverified spelling aliases are used to change assignments.

Excluding Ephydridae makes the 2023 littoral red terrestrial mean zero, so the
original proposed mean of eight log contrasts is undefined. The initial script
stopped at its positive-mean assertion; no pooled result from that attempted
specification was produced. Do not remove the red treatment or add pseudocounts.

Instead define a **different**, explicitly exploratory estimand: calculate
arithmetic mean captures across the fixed treatment menu first, then calculate
the aquatic-versus-terrestrial log response contrast. For 2023, average eight lit
cell means equally (four colours by two placements), and average the two dark
placement means equally. For 2020 use its original one lit and one dark arm.
Obtain sampling variance from independent within-cell covariance, with fixed
cell weights and delta derivatives at the aggregate means. Retain all zero
counts and all assigned treatments. Between-colour variability is not treated
as sampling replication. This is a standardized experimental menu mixture,
not an equal average of log effects and not a standardized common lamp.

Report the source-defined and family-derived versions as two representations
of the same two programmes, never as independent datasets. The two-programme
REML and modified Knapp-Hartung specifications above remain unchanged.
