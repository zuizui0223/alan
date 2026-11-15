# ALAN B | New verified diel fruit-set evidence and source-scale uncertainty (2026-11-15)

## What actually ran

Original Blareau & Requier source, DOI [10.1111/jen.70019](https://doi.org/10.1111/jen.70019), raw Figshare DOI [10.6084/m9.figshare.29930195.v1](https://doi.org/10.6084/m9.figshare.29930195.v1); first published 2025-09-12, issue volume 2026.

Previously source-verified SHA256 fruit_set = c5c4a169b345cccafa6d555cfa1726ae5f6606bb83312c2602bee3626cc2dd0f (80 inflorescence rows); fruit_quality = b988588bf14b73ca2f619845e09eedc92ef464b77e64417431fbb2b725334fd5 (223 formed-fruit rows).

**New executed source sensitivity:** [Actions run #43697496190](https://github.com/zuizui0223/alan/actions/runs/43697496190) SUCCESS, 6/6 specific tests pass. Source analysis scripts/analysis/alan_P1_smallfruit_inflorescence_sensitivity_20261115.py, tests/test_alan_P1_smallfruit_inflorescence_sensitivity_20261115.py, workflow alan-P1-France-850-inflorescence-sensitivity-20261115.yml, artifact #14154504067, ZIP SHA256 79d0cfa5aa25923ee54677ee8d3e595c16d9ed6b2f33b578955501bcccb0bdd9.

All four original daytime visitor access arms in each taxon-year were kept: Closed, Day, Night and Day+Night. Total 80 inflorescences, 850 initially assigned flowers and 258 formed fruits; 36 source plants, two survey years and different sites. A *species × year* is a descriptive stratum, NOT a new independently manipulated ALAN experiment.

## The actual original counts and direct marginal night-access dependence

V = 1 - (Day-only success / initially assigned flowers)/(Day+Night success / initially assigned flowers). No imputation or truncation. Day-only and open arms were assigned to distinct inflorescences, not the same flowers. V is about nighttime access and is NOT the effect of ALAN or identified moths.

| Crop | Year | Day-only successful fruits/initial flowers | Day+Night successful fruits/initial flowers | V |
|---|---:|---|---|---:|
| Ribes nigrum (black currant) | 2020 | 5/33 | 14/44 | **+0.5238** |
| Ribes nigrum | 2024 | 5/59 | 15/43 | **+0.7571** |
| Rubus idaeus (raspberry) | 2020 | 17/17 | 11/11 | **0.0000** |
| Ribes rubrum (red currant) | 2020 | 29/54 | 24/31 | **+0.3063** |
| Ribes rubrum | 2024 | 25/83 | 19/50 | **+0.2074** |

Source names were checked against original publication's total flowers: BC black currant 375, RB raspberry 54, RD red currant 421. Note that raspberry V=0 is a *boundary sample estimate* from 2 day-only and 2 day+night inflorescences, not a precise demonstrated absence of nocturnal pollination; authors report seed-per-fruit benefits for these same crops.

## Quantified robustness with correct limitations

One ORIGINAL inflorescence omitted from one of the two arms at a time (LOIO). Results below reflect only those leave-one group deletions.

| Crop/year | V LOIO min-max | Inflorescences in Day / Day+Night | 20,000 bootstrap resamples, **inflorescence-level descriptive percentile only** |
|---|---|---|---|
| Black currant 2020 | 0.015–0.630 | 2 / 2 | Suppressed (under 3 inflorescences/arm) |
| Black currant 2024 | **0.650–0.888** | 8 / 5 | **0.236–1.000**, fraction positive 0.9908 |
| Raspberry 2020 | 0.000–0.000 | 2 / 2 | Suppressed; degenerate 0 can still be uninformative |
| Red currant 2020 | 0.233–0.409 | 2 / 2 | Suppressed |
| Red currant 2024 | 0.053–0.347 | 7 / 6 | **−0.367–0.497**, fraction positive 0.8195 |

**These intervals are NOT valid plant-level 95% CIs.** Authors randomized pollination-access treatment by whole inflorescence (racemes), with multiple inflorescences per plant; original Figshare CSV lacks the plant ID and different-year site-level randomization, so covariance within plant cannot be recovered. Resampling inflorescences as independent is a sensitivity check, not an inferential plant-cluster bootstrap. The whole original French crop study counts as ONE publication, not five experimental ALAN replications. Even the 2024 black currant very positive LOIO does not prove biological effect modification across other populations.

## The important stage/join guard

Fruit-quality table has 223 measured-fruit rows, but the fruit_set table has 258 successfully formed fruits. Worse, **9 of 20 taxon×year×treatment cells have discordant successful-fruit vs fruit-quality-row counts, including 2 where quality rows exceed fruit counts.** Therefore 223/258 is not a simple random missingness proportion and cannot be treated as coverage. No reconstructed total seeds per initial assigned flower can be claimed without verifying individual fruit/inflorescence linkages and how quality files were curated. Separate seed-per-formed-fruit outcomes remain descriptive and conditional on post-treatment fruit formation.

## How this advances A without inflating it

B now gains a new authentic independent original *publication*, with five descriptive taxon-year night-access dependence strata and true source row counts. This is **not** five additional ALAN experiments. Do not transport V across countries/years to Macgregor, Knop, Li or Lu by genus/family similarity and call it a new meta-regression observation. The actual ALAN×V coefficient, valid 95% interval, and island-specific effect remain unestimated.

A decisive next dataset still requires genuine illuminated versus dark reproductive response plus independent day-night exclusion in the same experimental population, full unit-level n and shared-control covariance; Lu 2026 JPE has the necessary factorial design but source cell data are not publicly recoverable via the audited standard publisher links.

### Secondary exploratory thought (not causal)

An independent-channels approximation for open fruit probability, 1-(1-Pday)(1-Pnight)/(1-Pclosed), shows departures whose directions differ by taxon/year. This is only a mathematical baseline and should not be interpreted as evolved compensation or verified behavioral synergy without within-plant paired randomization and actual insect observations.

