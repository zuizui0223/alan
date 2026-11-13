# ALAN: Real paired field day/night pollination reproductive outcome versus true historical nightlight — global macro analysis

**Original source replay SUCCESS, 2026-11-13.** Replaces the non-identifiable between-study GloPL Sup/Bagout comparison in PR #50 with truly paired **day-versus-night pollinator-exclusion field fruit/seed production** from original published experiments.

[GitHub Actions real source #43380089481](https://github.com/zuizui0223/alan/actions/runs/43380089481); [aggregate result artifact #14015499514](https://github.com/zuizui0223/alan/actions/runs/43380089481/artifacts/14015499514). The previous PRs #43–50 and main are unchanged. This PR stays Draft because this is observational light exposure, not actual manipulated ALAN.

## Scientific estimand, authentic source and precise population

Published author experimental data: Kendall & Nicholson, *Pollination across the diel cycle: a global meta-analysis*, Ecology Letters 2025, DOI 10.1111/ele.70036. Original source 1,652 treatment rows, author Git study blob SHA 39a82d966afcdde670603add1c881ace8791ffe9. Published final 135 studies/139 plant taxa not equal to original source 140 study IDs; this source explicitly reconstructs experimental records.

**Do not pool original study×plant×metric effects across different field years or sites.** A coarsened source grouping gives 231 day-night fruit/seed study×species×metric pairs but mixes multiple originally reported field sites and years in some papers. The original corrected scientific source unit is published study × plant × fruit/seed metric × original measurement unit × exact original lat/lon × reported start/end experimental year × month bounds, both real source day/night treatment arms recorded together.

The independently authenticated source preflight [#43379308366](https://github.com/zuizui0223/alan/actions/runs/43379308366) froze **301 valid paired site/year/season/metric reproductive outcomes**, **118 independent published studies**, 123 originally named plant species, 135 original source field locations and 163 study×site×experiment-start-year units. Field years 1993–2022. Original source metrics: 178 fruit set, 123 seed set; 270 original paired outcomes described their source plants as native; 197 originally printed study coordinates, 98 coordinates inferred from place name and 6 from maps. Original elevation linked for all 301 source units. The author research source leaves END year unreported for 284/301 groups and reports multiyear date ranges for 17; exposure uses source START year only, never the study publication date or future light.

Original field-pair source SHA256 **21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af**; original source-year geolocation and source cohort size checked before any satellite-pixel join.

Original nightlight: **Li–Cao corrected PCNL v3 Zenodo record 17013414** with 30 distinct original GeoTIFF years **1992–2021**. Six independent 5-year annual sampling shards all completed SUCCESS under original published annual MD5, and a collector verified original raster years exactly once. At each original author experimental location, take mean log1p original raster values from 1992 to **year_start−1**. **All 301/301** paired day/night source units have full source-authenticated prior satellite data. Original joint source final SHA256 0ea518558117e517a73345a93daee132f9ba7f4525296a744df131df805a359f.

**Important support:** 193/301 = **64.1%** field experimental pairs had PCNL pre-experiment mean exactly ZERO. Only **9/118 independent studies** sampled >=2 original historical PCNL values with span >.03, and only **8** also used >=2 independently geolocated study sites. Therefore the historically prior light effect in this field-meta dataset is **primarily BETWEEN original studies and geography** even though each response compares day and night WITHIN a study. No field light switch-on or lamp-manipulation effect is identified.

## Original authentic external-study reproductive outcome

Within the same original study plant, true field site, experimental start year, months and original measured fruit/seed units, compare actual night-pollination-only versus day-pollination-only treatment outcomes by night_share = observed nighttime reproduction/(nighttime + daytime reproduction). Treatment means are weighted by original observation sample sizes where available. A night share of 0.5 is equal numeric day/night production; this is a **descriptive primary outcome**, not a reproduction standardized effect size or evidence of a significant within-species difference. Study-equal-weighted mean source night_share = 0.5701, source-median 0.5830. No direct pollinator taxon was resolved from study exclusions; original field daytime/nighttime work may involve non-moth taxa.

Models fit WLS with **EQUAL TOTAL WEIGHT FOR EACH ORIGINAL STUDY**, source-study clustered nominal 95% uncertainty; independent study × plant 2-way clustering as sensitivity. Original covariates in main source model: real experiment START year, geographic latitude/latitude²/longitude/longitude², independent author-published elevation, actual fruit versus seed reproductive metric. Original 20-degree region FE posthoc secondary model; source native, author original text coordinates, study-year precision, metric and pollination-hours source-quality sensitivity. All were predeclared **before observing PCNL correlation signs**.

## Global real original field reproductive result

Coefficients below in **percentage points of within-experiment nighttime reproductive share per one mean log(1+prior PCNL)**. Every CI straddles ZERO.

| Original field experiment source model | Observations / studies | β pp | Nominal 95% study-cluster CI pp | p |
|---|---:|---:|---|---:|
| **Main study/year/geography/elevation/metric adjusted** | **301 / 118** | **+0.3068** | **[−2.2388,+2.8523]** | **.8133** |
| Unadjusted prior PCNL | 301 / 118 | −0.2730 | [−2.5712,+2.0252] | .8159 |
| Additional original 20° regional FE | 301 / 118 | +1.7046 | [−1.0213,+4.4304] | .2203 |
| Source-native reproductive field plants only | 270 / 105 | −0.3569 | [−3.2409,+2.5271] | .8084 |
| Original precise coordinates *provided in text* only | 197 / 80 | +0.9575 | [−2.2300,+4.1451] | .5560 |
| Study groups without explicit multiyear range | 284 / 109 | +0.2575 | [−2.3405,+2.8556] | .8460 |
| Original **fruit-set** outcome alone | 178 / 100 | +0.3994 | [−2.3780,+3.1769] | .7780 |
| Original **seed-set** outcome alone | 123 / 71 | +2.5769 | [−1.0485,+6.2023] | .1636 |
| Independently measured day/night exclusion-duration normalized reproductive effect | 171 / 82 | −0.4014 | [−3.7082,+2.9054] | .8119 |

Original study×species TWO-WAY clustered uncertainty for main original field estimate: β fraction +0.0030676, finite-117-df 95% CI **[−0.0235093,+0.0296445]** (percentage points **[−2.3509,+2.9645]**), p=.8196. Same result.

**Primary interpretation: no robust/identified global negative nightlight association with directly measured relative nighttime reproductive benefit.** This is substantially closer to real pollination function than previous GloBI moth-plant or flower-photography proxies, and separates within-study day/night treatment comparison from between-study nightlight exposure. It is not a proof that artificial night lighting is biologically harmless.

## Important scientific boundaries

1. **Exposure** is historical corrected satellite nighttime radiance, not actual at-flower illuminance. A huge source fraction is exactly zero before experimenting. Averaging PCNL 1992–field-start-year−1 cannot establish lighting changes induced night pollinator declines.
2. **Confounding:** original experimental site geographic setting, vegetation, land use, author pollination-exclusion design, floral traits and intrinsic plant differences vary nonrandomly across the 118 published source studies. Only 8 studies have genuine original site/light contrasts suitable for even attempting within-study exposure separation; a new study-site FE model would be weak/non-identifiable, regardless of 301 paired outcomes.
3. **Different outcome and exposure designs** between the verified GloPL Sup-only pollen-limitation source and this genuine day/night-exclusion source mean coefficients are NOT directly subtractable, cannot be a mediation analysis and do not measure the same local plant populations.
4. **Not moth-specific:** day/night treatments allow plants to be visited by any activity-time-qualified pollen vectors; the author's meta paper explicitly does not provide standardized identities/effectiveness for nocturnal moths.
5. **Not island-only:** studies span true worldwide field coordinates including mainland and islands. No global island-specific light-syndrome effect was identified from this meta source.
6. **Original meta manuscript** uses Hedges' g phylogenetic multilevel study effect, and reports that roughly 90% studied plants lacked resolved diel success differences. Our original normalized paired night/day share is a different descriptive measurement, and does not claim to reproduce the article's standardized phylogenetic results.
7. The original day/night experiment end year is unreported for most records. Previous-light sampling uses reported field start date as the conservative time bound; that does not improve source date precision. Original coordinate quality is heterogeneous.
8. Several related exploratory models and source eligibility filters were examined post-hoc along the wider ALAN program. No single nominal p from this family qualifies as discovery of ALAN-caused biological function change.

## Consequence for ALAN research

Across past analyses, a negative source effect of 713 old-GloBI moth visitor plant labels was source-sensitive; strict independently measured moth-LIKE white/tubular floral syndromes did not yield a robust negative geography/genus-adjusted association; and a binary external day/night source plant-community comparison was weight/definition dependent. **The new real paired fruit/seed experimental benefit test also has no robust historical nightlight association.** Strong causal mechanism cannot be recovered by merely photographing more plants. Future macro work needs more reproduced local light-manipulation or repeat light exposure changes paired with effective pollen deposition, seed set and independent selfing assays of the same wild plants across islands.

## Reproducibility

- [Original source field/year pairing #43379308366](https://github.com/zuizui0223/alan/actions/runs/43379308366).
- [All original 30 PCNL annual MD5 raster years and study-cluster model #43380089481](https://github.com/zuizui0223/alan/actions/runs/43380089481), [final aggregate artifact #14015499514](https://github.com/zuizui0223/alan/actions/runs/43380089481/artifacts/14015499514).
- New scripts: original_DvN_fieldsite_year_source, authentic_diel_siteyear_PCNL_6shards, original_DvN_reproductive_vs_PRIOR_PCNL. Earlier source-corrected PCNL raster rehydration code and MD5 manifest copied byte-for-byte into separate PR #51 from the successfully audited GloPL source branch.
- Only aggregate model source results and separate annual-shard checksum receipts are published in final artifact; full original author study field coordinates and taxon×site/time rows are not redistributed in the user-facing output.
