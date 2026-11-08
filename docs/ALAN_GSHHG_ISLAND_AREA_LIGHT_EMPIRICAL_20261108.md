# ALAN | Does GSHHG island area modify the observed nightlight–floral-whiteness slope?

**2026-11-08 original-photo empirical exploratory extension — successfully executed in ALAN.**

## What has actually been estimated

Using the exact previously frozen original MOIS photo/WCVP, 2019/2021/2023/2024/2025 same-year VIIRS and R2024B UC WorldPop, and SHA-pinned GSHHG 2.3.7 L1 landmass geometry, the source-photo floral ROI9 white mass was regressed on same-year log1p(VIIRS) **within islands**, asking if the light slope varies with original GSHHG island area.

The predictor is **log10 island area (km²)**; therefore the moderation coefficient refers to **10-fold larger island area**. Area uses EPSG:6933 GSHHG polygon area; it is **not distance from mainland, geographic isolation, port development or a remotely sensed pollination measure**. In the strict core, 101 islands have at least one eligible photo-cell-year, but only **44 islands** with at least two distinct original 0.05° photo sites and distinct original light readings qualify for the within-island spatial analysis. The same site recorded in different years is NOT counted as distinct spatial locations. The broad measured-original-photo cohort includes 222 island landmasses, but only **95** qualify for the analogous within-island analysis.

Model in both preexisting, overlapping observed-photo cohorts:

```
white_pp ~ log_light + log_light:log10_area_c
        + log_pop + bio5 + log_rain + elevation + log_attempts
        + C(land_id) + C(year)
```

- `log10_area_c` is centered at the strict cohort's **equal-island median log10 area = 3.2272707936997804**. No category cutoff was used to select an advantageous size class.
- Clustered errors are at the **GSHHG island** level, repeated observation years retain year fixed effects. Nominal 95% intervals are **exploratory**, not multiplicity-corrected or experimental causal identification.
- Outcome is a **photographic original floral-region nine-palette white mass**, equally averaged across actually photographed accepted taxa per source photo cell; not verified petal-only morphology or native flora abundance.

## New genuinely executed coefficients

| Actual source-qualified contrast | Strict 5-photo/3-species cohort | Broader 3-photo/2-species cohort |
| --- | ---: | ---: |
| Source photographic island cell-years (all islands) | 710 | 2,026 |
| Distinct islands with any source-qualified cell | 101 | 222 |
| Source cell-years **within** selected multi-site islands | **650** | **1,884** |
| Distinct multi-site islands | **44** | **95** |
| **Nightlight × log10(island area) moderation (pp per tenfold area)** | **+2.2087283** | **+0.9983139** |
| Island-cluster nominal 95% CI | **[+0.0065948, +4.4108618]** | **[−0.3903994, +2.3870273]** |
| Nominal p | **0.04935** | **0.15679** |
| Nightlight slope at strict-median island area (pp per log(1+VIIRS)) | **−5.3501558** | **−1.3047489** |
| Nominal 95% CI for this conditional light slope | **[−10.1546161, −0.5456956]** | **[−4.0221250, +1.4126273]** |

**Result:** both measured-photo specifications give a **positive area interaction point estimate**. In the stricter sample, the nominal interval is only just above zero; in the broader source cohort its interval includes zero and the estimate is less than half as large. A positive interaction means that the nightlight–photographed-whiteness slope is *less negative* on larger islands, not that larger islands necessarily have more white flowers. The cohorts overlap and are **not independent biological replications**.

This is a **new exploratory area-structured spatial pattern**, not evidence that ALAN causes selection against white petals, a robust universal small-island mechanism, or the presence of a particular nocturnal-pollinator guild. Within-island fixed effects control all stable island-level differences but do not identify nighttime light causal exposure; they do not by themselves fix observer recruitment and taxon composition. VIIRS upward radiance is also not illuminance at any flower.

## Source support by island and equal-island modelling (completed in the same project)

A denominator audit found that equal photographic-cell weight was strongly dominated by the largest sampled islands even though geographic coverage spanned smaller islands:

| Real measured support | Strict | Broad |
| --- | ---: | ---: |
| Smaller-area third: independent islands / photo cell-years | 15 / 41 | 32 / 121 |
| Middle third | 14 / 64 | 31 / 171 |
| Larger-area third | **15 / 545** | **32 / 1,592** |
| Portion of photo cell-years from larger-area third | **83.85%** | **84.50%** |
| Islands with source-measured median across-site loglight span ≥0.05 | 36 of 44 | 84 of 95 |

To estimate a different, island-weighted target, a second model uses **1/n actual photo-cell-years within each island** as each cell weight, so the total weight of each original island equals one. This does not fabricate any new photos, and the regression formula and island-clustered errors are otherwise identical.

| Area × lognightlight interaction, in percentage points per tenfold island area | Strict 44 islands / 650 cell-years | Broad 95 islands / 1,884 cell-years |
| --- | ---: | ---: |
| Original **cell-year-equal** fit | +2.20873 [**+0.00659,+4.41086**], p=.04935 | +0.99831 [−0.39040,+2.38703], p=.15679 |
| New **island-equal** fit | **+5.15162 [ +2.55884,+7.74441 ]**, p=.000240 | **+3.38466 [ +1.31941,+5.44991 ]**, p=.001582 |

The island-equal approach gives a positive interaction with nominal intervals above zero in both pre-existing measured-photo cohorts. It thus provides **stronger exploratory evidence that the photographed-nightlight slope varies with island area** than the photo-cell-equal broad fit. The photo-cohort and island-equal estimates are different observational estimands; they are not four independent confirmations, and the models were introduced after inspecting existing results. The small islands also have fewer photographed sites, so equal island weighting gives sparse islands greater relative influence. Island area may be correlated with human activity, shoreline geometry, geography and non-random volunteer recruitment; no independent causal selection or pollinator mechanism has been identified.

The new, verified ALAN run is [#42439789459](https://github.com/zuizui0223/alan/actions/runs/42439789459) (**SUCCESS**), with original aggregate JSON SHA256 `35b30fb94b3a770142ac98498f835c7cec9c27f2273af41a5a98057763275810` and archived artifact ID `13600898564`. A durable source-verified numerical excerpt including the original JSON digest, model coefficients, island tertile denominators, and interpretive limits is committed at [results/alan_island_area_light_20261108_source_verified_summary.json](../results/alan_island_area_light_20261108_source_verified_summary.json). The original archival run #42439477013 remains a traceable precursor.

## Relation to the existing separate global result

The preceding source-verified **island vs mainland** global five-year interaction remains **−2.4970446 pp**, nominal 95% CI **[−3.9485472, −1.0455420]**, over 7,982 source photo cell-years and 710 island cell-years. This new **within-island area moderator** is a distinct estimand, and its strict regression uses only 650 eligible island cell-years. Do not algebraically combine the two coefficients or imply that area explains a percentage of the original island/mainland interaction. The geographically co-supported island/mainland comparison in the original work remains statistically unresolved.

## Execution and data lineage

- [Verified ALAN GitHub Actions run #42439477013](https://github.com/zuizui0223/alan/actions/runs/42439477013): **SUCCESS**, source inputs SHA-checked, original and synthetic tests passed, results artifact `alan-5year-GSHHG-island-area-nightlight-originalROI` (artifact ID **13601335880**; archived aggregate result JSON).
- Analysis: `scripts/analysis/alan_island_area_light_20261108.py`.
- Tests: `tests/test_alan_island_area_light_20261108.py`, including distinct-site vs repeat-year separation, fail-closed area classification/NA checks and known synthetic moderation.
- Workflow: `.github/workflows/alan-island-area-nightlight.yml`.
- Reuses old source data via historical signed MOIS Actions artifacts. No additional image downloads/inference; no unknown photo hue imputation, fake plant presence/absence or rejected missingness scenario.
- The output artifact has finite GitHub Actions retention; source archive provenance remains [ALAN migration record](../MIGRATION.md).

**Next biological test, not yet executed:** quantify **true mainland isolation** as a separate geographic variable (coastline separation, not merely island area). Do not silently treat small area as high isolation. The area moderator alone does not establish the channel-gated pollination hypothesis.
