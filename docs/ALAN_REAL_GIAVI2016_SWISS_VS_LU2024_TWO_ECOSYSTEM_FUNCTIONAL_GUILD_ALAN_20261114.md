# ALAN — independent Swiss field experiment, daytime channel turnover and cross-region transportability

**Original-source results reproduced: 2026-11-14.** Giavi, Fontaine & Knop (2021), Nature Communications 12:1690, DOI https://doi.org/10.1038/s41467-021-22011-8. Public original Swiss 2016 data: https://zenodo.org/records/4540407 . Original research archive includes data.Rdata, script.R and readme, all authenticated against original Zenodo MD5s. The original DAT frame contains **1,884 documented rows, 2,384 actually counted daytime visitor interactions, 21 plants, three orders, 37 paired survey occasions and six independent lit-vs-dark field pairs**.

[Original source manifest CI SUCCESS #43567987250](https://github.com/zuizui0223/alan/actions/runs/43567987250) and [completed raw field-guild analysis CI SUCCESS #43568387728](https://github.com/zuizui0223/alan/actions/runs/43568387728), aggregate artifact #14102984245. Data RData MD5 f60b21d82c163cbb53520fda0bc61985, SHA256 fb72e4e8a238371acacfb83e950348d9765eb24f72923802084c8fe5cff7e4bf.

## Independent Swiss 2016 vs Tibetan Lu 2024 experiment

This is a *cross-study sign and heterogeneity comparison*, NOT a pooled meta-analysis, as light, ecosystem, site design and original plant/pollinator taxa differ. Swiss was **daytime only**; Tibetan original study had day and night.

| Daytime guild | Tibet Lu 2026 actual lit/dark (2024 field) | Swiss Giavi 2021 actual lit/dark (2016 field) |
| --- | ---: | ---: |
| Diptera | +46.3% (527 to 771 visits) | **−3.339%** (569 to 550) |
| Hymenoptera | −20.0% (730 to 584) | **+11.159%** (466 to 518) |
| Coleoptera | −61.68% (167 to 64) | **+17.829%** (129 to 152) |

**Scientific reading:** ALAN response is not a fixed sign attached to an insect order. Giavi et al. had already published the significant plant × guild × light interaction; our incremental analysis is independent source-authenticated cross-study sign and within-site sensitivity, not a claim to discover the original result. Swiss plant abundance is recorded as source covariate plab but is not a synchronized census of available nectar-bearing flowers or direct pollen delivery.

## Original Swiss Diptera aggregate direction flips relative to typical field pair

Real six-field source aggregate **dark 569, lit 550**, net −19 (−3.34%), yet at **five of six** corresponding paired fields the lit total is HIGHER. Six anonymized site-pair original Diptera count directions: **+78, +23, +2, −149, +8, +19**. One extreme negative field determines the pooled source sign. The complete original leave-one-site-out check (all sites omitted in turn, NOT a selected omission) yields five negative pooled percent changes and **+34.30%** when that site is omitted. The complete-six-field median lit/dark ratio with source 0.5 continuity correction is **1.1751**, opposite the all-visit pooled ratio of **0.9666**.

Sampling uncertainty is substantial: **5,000 whole-field-pair resample descriptive 95% change CI [−48.42%, +60.19%]** for Swiss Diptera, so no stable cross-field causal direction. Corresponding Swiss Hymenoptera (+11.16%) whole-field resample [−45.87%, +74.51%], Coleoptera (+17.83%) [−35.45%, +115.63%]. No false independent-pair inference from the 1,884 raw plant-observation rows.

Source data selection: 21 botanical taxa observed, 15 with matched positive original visit evidence from both lit/dark arms in at least two separate field pairs and eight in at least five pairs. A missing plant-event row is unknown, not verified plant absence. No treatment effects on actual seed/fruit set are measured in this original source.

## Ecological relevance to syndrome/predator/island/mismatch programme

- Actual night fly rise in Tibet is not guaranteed to occur in another ecosystem or replace effective pollen delivery. Across 5 Tibetan independent blocks night Lepidoptera decline and fly increase were concordant, but the plant×block moth-majority designation failed transfer to heldout blocks (only 2/4 consistent).
- The Switzerland (2016) and Tibet (2024) experiments are independent primary research sites, but two source papers are not sufficient to estimate worldwide ALAN × floral syndrome × island effect.
- Important third evidence tier: preprint https://doi.org/10.64898/2026.07.10.737511 by Li et al. (2026), a UK **three-block light-configuration manipulation**, reports improved seed production at 25m array lights in three moth-associated phytometer plants, rather than uniform losses. The authors themselves note spatial arrangement differs from total emitted light and no moth movement was measured; treat light regime geometry as explicit hypothesis, NOT established corridor mediation.
- Giavi 2016 records no same-night flowering onset or effective pollen; it does not identify ALAN-induced seasonal/diel flower–pollinator phase mismatch, predator capture or predator floral pollen, or a verified island shift. The separate Swiss Grognuz experimental flowering raw database remained inaccessible by its Mendeley metadata API (HTTP 401), so no invented mismatch value.
- Night-vs-day functional reliance and true independently measured selfing ability must be acquired at local plant population level. True same-population island/mainland comparisons and predator two-role trials remain missing.

## Provenance and reproducibility

All three Zenodo files were independently retrieved and MD5 verified; the original study contains 1,884 repeated plant × site × visit-event RData rows, not 1,884 geographic replicates. Only six independent original light/dark meadow pairs. Analysis scripts:
- scripts/analysis/alan_Giavi2021_Swiss_original_data_source_audit_20261114.py
- scripts/analysis/alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114.py
- tests/test_alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114.py
- .github/workflows/alan-Swiss2016-original-sixsite-guild-reanalysis-20261114.yml

GitHub Actions #43568387728 success, aggregate-only artifact #14102984245. No plant×site×event original rows or site coordinates published, no redefinition of previous user results or main branch.
