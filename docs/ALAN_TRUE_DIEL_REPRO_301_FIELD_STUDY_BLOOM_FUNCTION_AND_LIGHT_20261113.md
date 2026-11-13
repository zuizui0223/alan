# ALAN: actual paired global day/night reproductive function, flowering time and historical nightlight

**Real original 301 field-study pairs, genuine author-published botanical bloom traits, original PCNL source and full stress tests verified 2026-11-13.**

Live original authenticated source & provenance: [Kendall & Nicholson (2025), Ecology Letters, DOI 10.1111/ele.70036](https://doi.org/10.1111/ele.70036); original authors Git study source blob 39a82d966afcdde670603add1c881ace8791ffe9 and trait blob 18381e550a7e405fe4f198909f0bb381bbecf1e0. Prior Li–Cao true source corrected PCNL from Zenodo 10.5281/zenodo.17013414, complete 30 MD5-checked actual raster years 1992–2021 at 1993–2022 original study field experiments, using only calendar years prior to original year_start. This is not local ground illuminance or manipulated artificial night lighting.

[Original field source and annual exposure replay #43380089481](https://github.com/zuizui0223/alan/actions/runs/43380089481); [source stress and botanical flower phenotype run #43387313576](https://github.com/zuizui0223/alan/actions/runs/43387313576). Aggregate-only artifact from the latter; no endangered study species×site coordinates or taxon identity rows published.

## Actual data: corrected independent units

Original 1,652 treatment observations are grouped by the **same published source study × plant species × reproductive outcome × original measurement units × exact original field site × actual experiment year and months**, then day-only and night-only pollinator-exclusion outcomes are compared as actual night reproduction / (night+day). The full source is **301 compatible site×year×outcome pairs / 118 independently published studies / 123 originally named plant species** from 135 original coordinates. 193/301 (64.1%) have zero preceding mean PCNL as measured by annual corrected satellite pixels. The original journal's final meta-analysis 139 species and Hedges g estimates should not be silently substituted for this reconstructive paired night-share descriptive outcome.

A conservative alternative collapsing within each study×species×outcome, requiring a unique original study field year and single source location, yielded 146 source comparisons / 87 studies /94 taxa, all 30/30 real-year PCNL MD5 sampled. That robustness source original β = **+0.984 percentage points** of night-share per mean prior logPCNL, 95%CI [−2.193,+4.162], p=.544, consistent with the larger source null and no evidence of negative direction. Its PR #52 was closed as duplicate/superseded by scientifically more complete PR #51; the rigorously reproducible numerical sensitivity remains documented in source Actions #43385562139.

## A. Historical artificial light and independently measured nocturnal reproductive function

Original study-equal weighted light association (adjust source experiment year, lat/lon polynomials, recorded elevation, reproductive metric; 118 independent original-study clusters): **β+0.3068 percentage points night-share per prior mean log1p PCNL**, CI **[−2.2388,+2.8523]**, p=.8133; two-way study+species clustered p=.820. This is a cross-study OBSERVATIONAL association between remote nightlight and experimentally observed day/night reproduction WITHIN each study.

### Source study and taxon leverage

All **118/118** one-publication deletions refit the original source. ALAN slope range in fraction: **−0.000859 to +0.007010** (−0.086 to +0.701 percentage points). **0/118 independent publication-omission CIs exclude 0**, smallest nominal p=.576. All **123/123** source taxon deletions: range **−0.001505 to +0.009361** (−0.150 to +0.936 percentage points), **0/123** CIs exclude zero; smallest nominal p=.506. Repeated leave-one-out estimates share nearly all the original source and are NOT 241 independent experimental replications.

### Actual historical illumination common-support limitations

| True original source circumstance | Observed support |
|---|---:|
| All original studies, paired site×year×metric outcomes | 118 studies / 301 source pairs |
| Zero prior PCNL | 193 source pairs from 69 studies |
| Positive prior PCNL | 108 source pairs from 56 studies |
| 10° regions in original field source | 53 |
| 10° regions with BOTH zero and positive previous satellite PCNL | 18 |
| **10° regions with ≥2 INDEPENDENT studies on EACH side** | **2** |
| **10° regions with ≥4 independent studies on EACH side** | **0** |
| 20° regions with both satellite source categories | 14/32 |
| 20° regions with ≥4 independent studies on BOTH sides | **3/32** |
| Original studies with real prior PCNL source variation >0.03 | **9/118** |
| Original studies with both actual dark and lit source cells | **7/118** |

This is the material reason even a study-equal source meta-regression **cannot infer a causal illumination effect**. Although the within-study night/day treatments are genuinely paired, PCNL contrasts are predominantly BETWEEN different studies/plants/geographies. Source black pixels do not prove darkness at the flower.

## B. Positive macroecological evidence: independent floral timing and real reproductive contribution

Original author botanical trait file independently labels original plant blooming periods as night, day, both or unclassified. Only source plants with **night or day documented** are used in the contrast:

- **196 actual field paired reproductive outcome site×years / 73 independent published studies / 74 plant species**.
- Documented NIGHT blooming: 140 paired reproductive observations, 48 studies, 46 species.
- Documented DAY blooming: 56 paired reproductive observations, 29 studies, 28 species.
- 28 botanical source families, **8 families with both classes**.
- Crucial: **only 4 independent published studies contain both night- and day-blooming plant species** (13 field response rows). This trait association is mainly BETWEEN botanical species and studies, despite within-study pairing of day/night pollinator exclusion arms. DO NOT interpret as experimentally induced bloom-time causation.

Night share = nighttime-only reproductive seed or fruit output / (nighttime-only + daytime-only output) within an actual published experiment, not an independent nocturnal pollen vector or effective moth-specific service.

| Source-verified night-bloom minus day-bloom difference in night share | Difference in percentage points | Study-cluster 95%CI (pp) | p |
|---|---:|---|---:|
| Original year+geography+elevation+fruit-vs-seed-adjusted | **+21.475** | **[+10.751,+32.199]** | **.0000867** |
| Add original 20° geographic region FE | +25.939 | [+14.846,+37.031] | .00000458 |
| Only original source-native field plants | +19.931 | [+8.519,+31.343] | .000619 |
| Only true site coordinates explicitly provided in paper | +28.121 | [+14.275,+41.968] | .0000688 |
| Original fruit-set outcomes only | +22.009 | [+9.541,+34.476] | .000540 |
| Site pairs with original day/night treatment hours available | +29.562 | [+17.626,+41.497] | .00000121 |
| Original botanical FAMILY fixed effects (8 mixed bloom families) | **+25.984** | **[+8.012,+43.956]** | **.00460** |

Study AND plant two-way cluster sensitivity on the original main β+21.475 pp: **CI [+9.411,+33.539] pp**, p=.000686; family-FE two-way cluster p=.00649. Therefore the positive direction survives distinct source methods of confounder control and study/species dependence, though not causal.

Original seed-set-alone lacks the predeclared independent study-count support and is honestly NOT ESTIMABLE.

**All 73/73 independent-study omissions** keep the night-bloom contrast positive: β range **+18.894 to +23.445 pp**. **All 74/74 single-species omissions** keep it positive: β **+18.894 to +25.409 pp**. Overlapping source deletions do not create new replication. More than one published plant and study contributes, but the limited within-study cross-bloom pairing remains a major confounder.

### Does nightlight change this trait-to-function relationship?

Interaction of original source blooming class × genuine previous-field-year PCNL in the same 196 paired reproductive outcomes: **β +0.017378 fraction (+1.738 pp) per prior logPCNL**, 95% study-cluster CI **[−0.049127,+0.083883] fraction (−4.913 to +8.388 pp)**, p=.6085; independent study+species two-way p=.582. **No supported light-specific decline or amplification of this real blooming-time/reproduction matching signal**.

## Scientific inference hierarchy

- **Supported descriptive macroecology:** a plant's independently documented nighttime blooming is associated with a substantially greater numerical nighttime contribution to fruit/seed production across 196 genuine paired field outcomes, and the positive direction is robust to source region, family, native-site labels, publication/species deletion and study/species clustering.
- **Not supported:** evidence that historical ALAN destroys, enhances or changes the actual local effective nighttime-pollination advantage. The field night/day comparisons are within an original study, but ALAN is mainly correlated ACROSS independent plants/studies, with only nine original studies with light variation and four studies comparing both day/night-blooming classes.
- **Not identifiable from this source:** actual moth-family pollinator identity, effective GloBI nocturnal service, local pollen deposition or stigma receipt, a selfing-mediated compensation response, island-specific native/wild plant decline, or evolutionary day-to-night transitions. The original bloom time may also be drawn partly from studies providing the outcome, so its ecological meaning must not be interpreted as a randomized contrast.

## Reproducible source and immutable evidence

Draft PR #53: new original-support and bloom-stability code:
- scripts/analysis/alan_DvN_real_field_light_study_influence_bloom_20261113.py
- scripts/analysis/alan_DvN_original_bloom_independent_robustness_20261113.py
- synthetic source tests and original SHA-rejoining Actions .github/workflows/alan-DvN301-original-study-leverage-bloom-20261113.yml.

The CI reuses original 301 authorized source site/year original artifact from run #43379308366, and all six original MD5 raster shard results from #43380089481 instead of re-downloading 30 year rasters; the source 301 field study paired outcome and year exposure are independently rechecked before any new p-value. Final output is aggregated-only. Earlier PR #43–51 and main remain intact.
