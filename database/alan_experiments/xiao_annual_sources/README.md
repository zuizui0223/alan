# Annual herb: reproductive timing and reproductive proxies

Reviewed source: [bioRxiv v2](https://doi.org/10.1101/2022.12.11.519667),
posted 2025-08-21. The bioRxiv API links it to
[the final paper](https://doi.org/10.1016/j.biocon.2025.111472).
The final paper's complete Methods/Results have not been verified.
These versions are one publication family, not independent experiments.
The indexed version-link table currently requires both members in the search
export; this supplementary link is recorded separately in `version_link.json`.

The field study used Elsholtzia densa in 20 open plots, 10 per light arm,
within a larger 40-plot experiment. Individual plants are subsamples.
Source coordinates: 50.17 N, 119.39 E; 2021; 6500 K LED; reported mean
21.15 lux. These coordinates and a grassland description do not establish
site wetness. This is distinct from the greenhouse Xiao herbivory study.

## Recorded numerical evidence

`reported_model_coefficients.csv` transcribes **all 14** ALAN coefficients
from v2 Table 2, PDF p.32, retaining both 95% and 90% Bayesian credible
intervals. The rendered table was visually inspected. The four phenophase
coefficients and fruiting-inflorescence coefficient received an independent
numerical check in the [design review](../search/xiao_annual_herb_design_review.md).

This is a supplement of author-reported model results, not a reproduction
of their models, raw arm summaries, or new sampling-variance lnRRs. No SE
was inferred from Bayesian intervals. Non-Gaussian link functions are not
explicitly specified in the reviewed Methods; the coefficients remain on
their reported scales without exponentiation. The first trait is square-root
transformed. Coefficients, prose percentages and arithmetic group differences
must not be substituted for one another. In particular, blooming's model
coefficient is -5.494 days, whereas the prose gives 6.2 days of advancement.

Reproductive traits include seeded-inflorescence counts and undeveloped top
inflorescences, not direct seed counts or fertilization success. The fruiting
inflorescence coefficient's 95% interval crosses zero, although its 90%
interval does not. Pollinators, pollen delivery and phenological mismatch
were not measured. Flowering onset was defined by the first purple
inflorescence; preserve this definition when combining phenology studies.

## Remaining extraction requirements

The preprint states that data are available on request; no public raw table
was recovered. No author contact was made. Final-version verification,
plot-level arm means and dispersion, treatment allocation, and cross-trait
covariance remain necessary before admission to the planned raw-data meta-analysis.
Do not count the 14 traits or 200 marked plants as independent studies.
No main-catalogue rows were added.

`source_manifest.json` records hashes of recovered PDF, XML and version API
metadata. External cache: `work/gbif_annual/alan_experiment_sources/xiao_annual/`.
