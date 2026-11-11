# GloBI all-positive records × 2019–2025 island/mainland ALAN: equal-month source audit

**Actual original-source observational analysis executed on 2026-11-11. Not a causal plant loss or local moth-pollination finding.** Neither fd nor island Chapter 1 was modified.

## Why the annual reconstruction was required

The earlier five-year 2019/2021/2023/2024/2025 source used monthly VIIRS from January, April, July and October (at least three available). In 2022, all 52,075 source points have an October provider error (PROVIDER_RASTER_400), not zero radiance. In 2020, WorldPop Global1 is not comparable to the R2024B unconstrained source elsewhere.

We recomputed all seven years using **the same January/April/July source months only**, requiring at least two valid monthly values. We retained provider-OK finite negative background-corrected source radiance as valid, clipping values only before the log1p brightness transform; the 3-of-3-month gate is a separate sensitivity. The 2020 population product was not substituted: the new model **omits population for ALL years**. This is still three sampled months, NOT an official annual VNL composite.

## Authenticated source and denominator

- Original one-photo/observation canonical SHA256: f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c
- Original WCVP exact accepted-ID mapping SHA256: dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2
- WCVP accepted names SHA256: 9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3
- Historical bounded GloBI positive moth/butterfly interaction source: ALL **61,967** existing rows, SHA256 a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994. This is NOT the full GloBI v0.11 2.5 GB dataset.
- Original source geography/monthly VIIRS SHA256 receipts for every year are recorded in the numerical source audit and bundled reproducibility code.

There are **284,950** photo-site-years with an accepted WCVP taxon in source scope. **48,005** pass >=3 photos, >=2 actually photographed plant species and at least two of Jan/Apr/Jul VIIRS months. The 2020 environment source is **only the archived final 30,630 of 37,084 original known-land 2020 cells**, not a complete 2020 census.

2020's direct same-year matched island/mainland sample: **83 island / 505 mainland** photo-site-years in 23 independent geographic 5-degree blocks. 2022: **181 island / 1,876 mainland** in 35 blocks. Earlier quick count could differ by a boundary site because this reconstruction consistently uses the original grid encoding floor((longitude+180)*20). The earlier premature decision to exclude negative but provider-OK radiance values was corrected before these final 7-year models.

## Model and real results

Outcome: percent of distinct actually photographed accepted plant species in a photo-cell×year that have a **positive GloBI moth-family flower-interaction record somewhere**. An association elsewhere is NOT an observed local moth visitation or effective nocturnal pollen transfer.

Compare only shared observation-year×5-degree geographic groups having BOTH island and mainland photos. Balance island/mainland total site weights within those groups, adjust log photo attempts, log documented photographed species count, BIO5, rainfall and elevation; geographic 5-degree blocks clustered across years. This is a new estimator and a new cohort versus the previous population-adjusted four-month model.

| Observed source scope | Island site-years / total compared | 5-degree blocks | Difference in island−mainland nightlight slopes (pp/loglight) | nominal 95% CI | nominal p |
|---|---:|---:|---:|---|---:|
| Previous 5 observed years recomputed on identical 3 source months | 950 / 10,989 | 79 | **−1.877** | [−3.681, −0.073] | .0416 |
| Six years, with 2022 but exclude incomplete 2020 | 1,131 / 13,046 | 80 | **−1.667** | [−3.351, +0.017] | .0523 |
| Seven years with **partial 2020** | 1,214 / 13,634 | 82 | **−1.203** | [−2.843, +0.437] | .1482 |
| Seven years, Jan/Apr/Jul ALL valid | 1,015 / 12,115 | 70 | **−0.682** | [−2.383, +1.020] | .4269 |
| Seven years, at least 5 photos AND 3 species | 352 / 3,496 | 41 | **−0.153** | [−2.273, +1.967] | .8849 |

The 7-year linear calendar change in the slope difference is **−0.448 pp/year [−1.230,+0.334], p=.258**, not statistically resolved. The negative association is not robust to including the additional years or strictly matching photo/month coverage.

Repeated-site panel: **2,238 authentic photographic sites, 6,453 site-years, 556 island site-years and 34 geographic blocks**. Site FE slope difference **−5.018 [−24.711,+14.676], p=.608**; calendar change **−0.275 [−1.159,+0.608], p=.530**. There is no reliable time trend.

Auxiliary observed photo-species moth-minus-butterfly record contrast: year-dependent slope difference **−1.083 pp/year**, nominal p=.00079, but iNaturalist GloBI source exclusion gives **−0.400 [−0.944,+0.144], p=.147**. This is source-sensitive evidence, NOT observed biological moth decline.

## Critical boundaries and next gate

The 2020 partial source is not assumed spatially representative. Species missing from photographs or GloBI are not proven biologically absent. Present-day GloBI positive interaction records do not demonstrate local moth dependence, reduction in pollination, naturalization failure or ALAN causation. All tests and associations here are exploratory.

Restore the earlier 2020 known-land numeric partitions. For a defensible longer environmental trend, obtain genuinely **annual** EOG VNL from the official annual VIIRS sources, audit cloud-free coverage (cf_cvg) and distinguish V2.1 earlier-year from V2.2 2022–2024 source processing. Official catalog: https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_ANNUAL_V22 . The reproducibility ZIP supplied with this analysis contains the independently tested full source builder, model script, numerical aggregate tables and 7 passed synthetic mathematical tests. The original species×location data are held locally and never committed.
