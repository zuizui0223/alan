# ALAN A | Knop 2017 original factorial fruit experiment: true source grain (2026-11-15)

**Original paper:** Knop et al. (2017), Nature 548:206–209, DOI https://doi.org/10.1038/nature23288 . Methods and figure verified in the author-uploaded complete manuscript (https://www.researchgate.net/publication/318862971_Artificial_light_at_night_as_a_new_threat_to_pollination). Important: the independently light-assigned field units are meadows, not the numbers of capitula or individual achenes.

## New certified denominator audit

| Component | Actual source design | Relevance |
|---|---|---|
| Nocturnal flower visitation network | 14 meadows, **7 illuminated + 7 dark** (2015); additional pre-light day/night observations 2014 | Not 14 reproductive experiment clusters |
| *Cirsium oleraceum* fecundity experiment | **5 illuminated + 5 dark = 10 meadows** in 2015, with 5 matched site pairs | Only **10 distinct light-assigned field units** for fecundity |
| Plants | 10 plants per meadow = 100 total | 5 nearby plant pairs per site, one randomly bagged in each pair |
| Factorial design | dark/open, dark/bagged, lit/open, lit/bagged | **25 plants per each of 4 combinations** |
| Plant reproductive response | TWO capitula per plant, counted developed white fruits versus aborted brown ovules | **50 capitula per combination**, 200 total; capitula are nested in plant and meadow |
| Source models | Poisson developed fruit count and binomial developed-fruit proportion per flower head, cluster random effects and overdispersion | Response is **achene number/proportion per capitulum**, not binary successful fruit per independent plant |

The original figure states n=50 capitula per factorial cell, but this is NOT n=50 independent plant units per cell and certainly not 50 independently illuminated meadows. Raw original flowering observations, site-level effects, and exact variance components are **not openly downloadable from this source**: original Methods state data and code available from corresponding author on request.

## Published original patterns and what they do not prove

On open-to-visitors plants, developed fruit count decreased under ALAN (author-reported z=-2.5, p=.014). Bagged plants' light difference was unresolved (author z=+1.47, p=.140); open fruits fraction has author p=.018 in Fig 4b. Reported p values are original model tests, **NOT meta lnRR standard errors**.

Crucially **one treatment arm showing p<.05 and another p>.05 does NOT establish a significant light×bagging interaction**. The original model considered interactions with bagging and site/pair effects and selected terms by AIC. Reconstructing its treatment-by-bagging coefficient and covariance requires the Extended Data Table 3 numerical rows or plant/site data; a new indirect pollinator mediation coefficient must not be fabricated.

Bagging suppresses *all* animal flower visitors (not selectively nocturnal insects). Bagging vs open cannot supply night-only marginal dependence V = 1 - reproductive(dark,day-only)/reproductive(dark,open). Original day/night flower visits cannot replace a reproductive exclusion experiment. Thus this is an important **direct experimental ALAN fruit response** but is **not by itself a paired observation for the core night-dependence meta-regression**.

Site structure matters: the original fruit experiment has 5 illuminated/dark site pairs, while 2015 nocturnal visits drew from 7 site pairs and 2014 daytime/nighttime visitor comparisons are earlier than the light experiment. Avoid cross-year or extra-site joins presented as same plant reproductive dependence.

## Inference / priority

The original 13% fruit-set reduction is an authors' summary from this complex nested design; do not treat it as a precise independent-flowers logRR with n=200. The full data have not been recovered, so this is a **design/grain certified, effect variance not certified** study.

Prioritised data request: matched site-pair ID, site treatment, bagging/open, individual plant ID, head ID, number developed fruits, total viable fruit plus aborted ovules, plant biomass, herbivory, observed Lux, statistical model output from Extended Data Table 3 (particularly treatment-by-bagging β/SE if retained). Aggregates at true 10-meadow level are an acceptable alternative to exact plots/coordinates.

Keep this source as ONE original reproductive publication, not 14 experiments, 100 repetitions, 200 independent treatments, or a positive nighttime-dependence slope test.

Repository: config/alan_A_reproductive_screen_20261115.json source_audit.knop2017 (added); tests/test_alan_A_knop_design_grain_20261115.py; GitHub Actions validation below.
