# Shared-light-group descriptive outcomes, 2026-11-20

Source: [Dananay & Benard 2018 original workbook](https://doi.org/10.5061/dryad.pq54gh5), checksum-verified Zenodo mirror4932781; methods evidence and unresolved eligibility are in [review](../search/amphibian2018_review.md). This is a descriptive extraction, not a new meta-analysis or adopted effect catalogue entry.

## Estimand fixed for this extraction

Each source-labelled spatial group contains four mesocosms with the complete lid-by-toad treatment combination. Average these four equally, then average the four illuminated and five natural-light groups separately. This is an average across that experimental treatment mixture, not an unrestricted wetland-community effect. The supplementary design states lamps were centred within illuminated groups; workbook Metadata defines Block as a spatial grouping of four mesocosms. This supports group-level descriptive aggregation, but no mesocosm-ID physical lamp diagram was recovered. Original subgroup values remain in the archived workbook.

Use the last recorded algal date (Leaf2July, Water2July), not an average of June and July. PredatorTotal and NonPredatorTotal are totals across both pipe-sampling days according to Metadata; they are not final-time snapshots. Retain endpoint names: the NonPredatorTotal metadata confusingly says 'non-amphibian predators', so no independently validated trophic classification is claimed.

|Source endpoint|Illuminated group mean|Natural-light group mean|Descriptive ratio|
|---|---:|---:|---:|
|Leaf2July|49.1875|79.79775|0.616402|
|Water2July|21.6158125|20.9569|1.031441|
|PredatorTotal|5.6875|4.55|1.25|
|NonPredatorTotal|1.5|0.85|1.764706|

Algal endpoints are the workbook chlorophyll/relative fluorescence indicators, not dry biomass or plant reproductive performance. Invertebrate counts describe colonization/sampled presence, not attacks, consumption or predation rate. Different response dates and sampling windows preclude interpreting this table as a contemporaneous mediated pathway. Lids alter light as well as colonization. Do not infer causal mediation from opposite directions of these marginal means.

## Verification and limits

`python database/alan_experiments/amphibian_sources/summarize_groups.py` writes36 group-endpoint rows and4 arm summaries. It verifies source MD5, complete one-to-one36-ID joins, light agreement, four factorial combinations per group and no missing endpoint values. No missing records are converted to zero. Group counts are4/5, not16/20 independent lights; responses from the same groups are dependent. The initial descriptive extraction made no interval or lnRR-variance estimate; the conditional extension below adds these without main-meta adoption or an independent-study count increment. Final source and indoor-phase review remain incomplete; an unsuccessful EuropePMC retrieval and an NCBI BioC HTTP200 error page did not supply primary text.

An independent reaggregation reproduced all four means and ratios, verified complete2x2 treatment combinations and no missing observations, and confirmed the trophic-label ambiguity. Workbook Metadata documents exclusion of mesocosms13–16 because of shade; they are not reconstructed as observed zeros. The original metadata gives no measurement units for the two chlorophyll columns, so none are assigned here.

## Conditional uncertainty extension

The [specification](UNCERTAINTY_PLAN.md) was written before interval calculation. Using four lit and five natural-light groups, delta-method sampling covariance and endpoint-wise Welch-t intervals gives:

|Endpoint|lnRR|Sampling variance|Ratio|Approximate95% interval|
|---|---:|---:|---:|---:|
|Leaf2July|-0.483856|0.060078|0.616402|0.340916–1.114501|
|Water2July|0.030957|0.037117|1.031441|0.651886–1.631990|
|PredatorTotal|0.223144|0.127460|1.25|0.528153–2.958424|
|NonPredatorTotal|0.567984|0.313197|1.764706|0.392342–7.937422|

Every interval includes1. This does not demonstrate equivalence or absence of effects. The estimates assume the source-labelled spatial groups are independent illumination units; a physical ID-to-lamp diagram has not been recovered. Four and five clusters and sparse counts limit the delta approximation. Full4x4 covariance is retained in `conditional_covariance.csv`; the four endpoints must not become four independent studies or be averaged into a common ecological response.

This terminal-date algal estimate differs from the published repeated-date model and its treatment of experimental units. It is neither an exact reproduction nor a refutation of that model. These conditional effects remain outside the main108-effect catalogue pending source/design eligibility adjudication. `compute_conditional_uncertainty.py` checks covariance symmetry/positive semidefiniteness, scalar variance agreement and preservation of all descriptive ratios.

Independent verification reproduced all four lnRRs, scalar variances, Welch degrees of freedom, intervals and all16 covariance entries. Scalar/matrix diagonal discrepancies were below7e-18; minimum eigenvalue0.015970. This validates arithmetic under the stated assumptions, not final adoption.
