# Capture-origin exploratory synthesis: independent audit

## Primary-source design

[Parkinson2020](https://doi.org/10.1371/journal.pone.0240138) used20 independently illuminated/unilluminated mesocosms at Haven Hill Lake; pan captures came from one August2016 night. Failed fish treatments were omitted by the authors. Both invertebrate origins were measured in each trap. [Parkinson2023](https://doi.org/10.1098/rstb.2022.0361) randomized30 stations at Waterfowl Pond on one June2022 night: four colours plus dark, three stations per shoreline placement. Violet intensity differed from other colours. Pan catches include spiders as well as insects. Neither experiment measured emergence production or population growth.

The figure legends assign Ephydridae to aquatic origin in2020 but terrestrial origin in2023. Author-defined origin categories therefore differ. Both experiments share investigators and region, but have different lakes, years and physical units. Independent sampling is plausible; broad geographic or research-group replication is absent.

## Independent data and estimand audit

Inspected original2020 pan CSV,2023 taxon/trap tables, extraction scripts, the9-row contrast table and full covariance. Exact lowercased source-spelling intersection contains15 families; Ephydridae is the only origin disagreement within that intersection. No spelling variants or higher/lower taxonomic ranks were silently reconciled.

Ephydridae counts:

| Programme/placement | Dark | Green | Violet/purple | Red | White/ALAN |
|---|---:|---:|---:|---:|---:|
|Parkinson2016|4|—|—|—|11|
|SevenPonds2022 littoral|0|2|2|3|3|
|SevenPonds2022 riparian|0|0|0|0|2|

The2020 family-column sums also differ from its supplied aquatic subtotal:

| Mesocosm | Treatment | Sum of aquatic-labelled families | Supplied aquatic total |
|---:|---|---:|---:|
|5|no ALAN|14|12|
|7|ALAN|170|172|
|8|no ALAN|19|20|
|18|ALAN|182|183|

All other aquatic unit sums match. Family-based aquatic total1882 versus supplied1884; terrestrial columns sum108 and match every supplied terrestrial subtotal. Preserve both source representations; these discrepancies do not justify guessed repairs.

## Pooling decision

A source-defined treatment-menu average is mathematically computable, but cannot be presented as a harmonized biological aquatic-versus-terrestrial effect. Before promoting a pooled estimate, prefer a clearly separate family-derived analysis that excludes the discordantly classified Ephydridae from both origin groups in both programmes. Sum the remaining source-labelled families directly, retain original source-defined estimates, and document unresolved spelling/rank distinctions. Do not merely subtract Ephydridae from inconsistent source totals and call the result fully reconciled.

Within-programme equal weighting defines a geometric mean of response-ratio ratios over the observed treatment menu. It is not the response ratio of pooled catches and not a common white-light exposure. SevenPonds averages8 dependent contrasts, not8 independent replications. Preserve shared controls and within-trap origin covariance through a'Va; never use sum(diagonal)/64 alone.

Independent checks on the current source-defined table give Parkinson2016 mean1.6305059783/variance0.1104314730 and SevenPonds2022 mean1.8686285252/variance0.1761166998 (ratios5.10646 and6.47940). These numbers describe the existing source-defined table, not the proposed harmonized variant.

Two-programme REML plus modified Hartung–Knapp max(1,q), df1 is permissible only as an explicitly provisional supplementary summary of these two menus. It does not estimate a wet/dry moderator, publication bias, causal trophic mediation or a transportable species-general effect. A zero heterogeneity estimate cannot establish homogeneity. With n3/station cells, sparse counts and possible shared attraction, delta-method variance and no-interference assumptions remain limitations; modified intervals do not repair measurement or classification problems. No model was fitted by this reviewer at this stage.


## Revised arithmetic-mixture estimand and numerical verification

After the taxonomy audit, excluding Ephydridae leaves all three littoral-red terrestrial catches zero. Therefore the original average of eight log contrasts is undefined for the restricted variant. Dropping red or adding pseudocounts would change the evidence. The amended estimator instead averages **counts first**, with fixed equal cell weights (eight lit cells and two dark habitat cells in2023), and then computes `log(A_lit/A_dark)−log(T_lit/T_dark)`. This is defensible as the relative response of the deliberately specified mixture of sampled lighting/placement conditions. It is not an equal average of log contrasts, a universal species response, a common exposure, or the response of a typical taxon. The amendment follows inspection of the data and must remain labelled exploratory.

Independently re-read the original2020 CSV and parsed all30 blocks of the original2023 Excel `raw site data` sheet with openpyxl. Origin assignments followed the published figure and explicitly identified source spellings. For each cell, computed the2×2 sample covariance of aquatic/terrestrial counts divided by its station count; summed cell matrices with squared fixed weights. Applied the four-variable log-ratio gradient to the resulting lit/dark mean covariance. This calculation did not reuse the implementation's centered-score formula or precomputed2023 trap totals. All four programme estimates and variances agreed:

| Variant | Programme | Log origin-response contrast | Sampling variance |
|---|---|---:|---:|
|Source-defined counts|Parkinson2016|1.6305059783|0.1104314730|
|Source-defined counts|SevenPonds2022|2.4664890643|0.3492948553|
|Family counts excluding Ephydridae|Parkinson2016|1.6341370798|0.1114596322|
|Family counts excluding Ephydridae|SevenPonds2022|2.4994590835|0.3505350917|

REML was independently optimized using the sole residual contrast likelihood, `log(v1+v2+2τ²)+(y1−y2)²/(v1+v2+2τ²)`. Modified Hartung–Knapp variance and t1 interval were then calculated separately. Numerical agreement with the saved outputs was confirmed to relative tolerance1e-7 (tight optimizer tolerance was required):

| Variant | REML τ² | Pooled response-ratio ratio | Modified HK95% interval |
|---|---:|---:|---|
|Source-defined counts|0.119570696|6.72366744|0.0457023–989.177846|
|Family counts excluding Ephydridae|0.143393723|6.88025420|0.0375877–1259.397280|

**Verification: passed for the revised calculation.** Both intervals span1 very widely. Positive programme estimates do not establish a precise population-average effect. The two variants are dependent analyses of the same samples, not independent replications. The classification audit, source subtotal conflicts, unequal treatment menus, possible trap interference and limited geographic/research-group replication remain substantive limitations.

Files verified: [implementation](../synthesize_capture_origins.py), [programme estimates](../capture_origin_synthesis/programme_estimates.csv), [pooled estimates](../capture_origin_synthesis/pooled_estimates.csv). Original records: [2020 author CSV](https://github.com/LizzParkinson/ALAN_Data/blob/6985359fd29dfe321c740691bf97ea9d564db769/Aquatic%20and%20Terrestrial%20Invertebrate%20Inputs.csv), [2023 original supplement](https://doi.org/10.6084/m9.figshare.c.6856475).
