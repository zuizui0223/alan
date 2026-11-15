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

**Correction after obtaining the original Extended Data Table 3:** An open-treatment p<.05 and bagged-treatment p>.05 do NOT, by themselves, establish a significant interaction. However, **the original authors DID directly test the lighting × bagging interaction.** In the indexed archived Nature Table 3, the developed-fruit-per-head fraction (binomial link) reports interaction **beta=-1.40, SE=0.592, p=0.018**; the developed-fruit count model reports interaction p=0.014. This is real source-reported evidence of differential **all-visitor access** responses under ALAN, not proof of a nocturnal-only mechanism. The binary percentage interaction has a ratio-of-odds-ratios exp(-1.40) about 0.247 if the extracted factor direction is as labelled; **not an ALAN risk ratio.** Its Wald-only 95% interval exp(-1.40 ± 1.96×0.592) is about 0.077–0.787, with source-link-scale model assumptions and field clustering as in the authors' fitted GLMM.

**Verification caveat:** the PDF as indexed at [public archived extended table](https://planesynormas.mma.gob.cl/archivos/2021/proyectos/ANTECEDENTES_ENTREGADOS_POR_ROC-OIKONOS.pdf) gives count-model interaction **+0.76, SE=0.309**, which is directionally inconsistent with the original Figure 4 when naively read as treatment Lamp minus Dark; the original table cell imaging or factor contrast coding must be inspected before signing this count estimate. The original archive is a reproduced PDF, not a separately deposited modelling dataset. Its 11.4MB downloaded copy was accessible in Actions but this page had no machine-extractable original text for independently confirming the sign. We **do not** average the count and fraction interaction coefficients or misrepresent them as paired day-night exclusion. Sample-level covariance and the actual model's count and fraction y/SE remain unavailable.

Bagging suppresses *all* animal flower visitors (not selectively nocturnal insects). Bagging vs open cannot supply night-only marginal dependence V = 1 - reproductive(dark,day-only)/reproductive(dark,open). Original day/night flower visits cannot replace a reproductive exclusion experiment. Thus this is an important **direct experimental ALAN fruit response** but is **not by itself a paired observation for the core night-dependence meta-regression**.

Site structure matters: the original fruit experiment has 5 illuminated/dark site pairs, while 2015 nocturnal visits drew from 7 site pairs and 2014 daytime/nighttime visitor comparisons are earlier than the light experiment. Avoid cross-year or extra-site joins presented as same plant reproductive dependence.

## Inference / priority

The original 13% fruit-set reduction is an authors' summary from this complex nested design; do not treat it as a precise independent-flowers logRR with n=200. The full data have not been recovered, so this is a **design/grain certified study with a published all-visitor lighting interaction β/SE**; nevertheless the desired open-lit/open-dark fruit log-RR and the independent night-dependence moderator V remain unverified.

Prioritised data request: matched site-pair ID, site treatment, bagging/open, individual plant ID, head ID, number developed fruits, total viable fruit plus aborted ovules, plant biomass, herbivory, observed Lux, original complete model interaction output and factor-level/contrast coding to independently confirm Extended Data Table 3, plus within-site covariance for the actual open light/dark contrast. Aggregates at true 10-meadow level are an acceptable alternative to exact plots/coordinates.

Keep this source as ONE original reproductive publication, not 14 experiments, 100 repetitions, 200 independent treatments, or a positive nighttime-dependence slope test.

Repository: config/alan_A_reproductive_screen_20261115.json source_audit.knop2017 (added); tests/test_alan_A_knop_design_grain_20261115.py; GitHub Actions validation below.
