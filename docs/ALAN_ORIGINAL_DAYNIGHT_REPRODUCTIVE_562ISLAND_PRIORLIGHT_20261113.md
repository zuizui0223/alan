# ALAN: Original daytime vs nighttime reproductive exclusion × prior nightlight across 562 islands

**Original-source analysis executed, 13 November 2026.** Source-locked full-photo flora and author treatment data; exploratory observational cross-study composition, NOT an effect of ALAN on local reproductive output.

Primary replay: [Actions #43342437842](https://github.com/zuizui0223/alan/actions/runs/43342437842) SUCCESS; extended predeclared source-sensitivity replay [Actions #43342907748](https://github.com/zuizui0223/alan/actions/runs/43342907748) SUCCESS; [result artifact 13999908320](https://github.com/zuizui0223/alan/actions/runs/43342907748/artifacts/13999908320).

## Original independent functional evidence

Author data: Kendall & Nicholson 2025, Ecology Letters, DOI 10.1111/ele.70036. Git-blob pins:
- day/night study raw source: 39a82d966afcdde670603add1c881ace8791ffe9
- species-level bloom/trait source: 18381e550a7e405fe4f198909f0bb381bbecf1e0.

Exactly 1,652 original treatments, 262 study×species×outcome records with both daytime and nighttime pollination treatments. Of these, 259 have strictly comparable finite positive combined responses and aligned units; **228** directly compare fruit set/seed set, spanning **131 independent source studies and 134 distinct raw botanical names**. Experimental reproductive records are NOT guaranteed effective moth-only pollination; night can include bats, birds, non-moth insects, or other pollen pathways.

Within matched study/species/reproductive outcome, create night_fraction = night_measure/(night_measure + day_measure), dropping incomplete/zero-sum endpoints and incompatible units. Mean across metrics within study/species; median across source studies within plant species. This is a descriptive source classification, NOT a standardized effect-size meta-analysis, not a p-value for individual plant species. In original study plants, 75 have greater numerical reproductive response under night-only treatment, 59 have greater daytime response. Distinct source treatment duration sensitivity uses time-normalized response where both day/night survey hours valid.

## True original island source and study overlap

Original 93,209,351 photo-ID index; accepted WCVP photo flora 14,330 species; 16 original source groups and 33-year harmonized 1992–photo-year-minus-one Li-Cao PCNL; **562 exact GSHHG islands / 77,647 0.05° cell×observation years / 2,012,925 unique species×site×year observations**.

Only **52** of 134 source experimental species have exact accepted scientific name matches in the original photographed island flora, **27 night-numerically-stronger** (recorded on 288 islands/10,608 site-years), **25 day-numerically-stronger** (285 islands/6,534 site-years). An unobserved functional class is UNKNOWN rather than day-pollinated, not obligate moth and not genuinely extinct.

The exact same photo-cell×year shows BOTH source documented classes on **2,565 cell-years across 130 islands**; the narrower within-island-year light-gradient eligible subsample is **2,373 cell-years, 39 independent islands and 111 island×years**.

## Main signed source coefficients (all observational, post-hoc)

Response: photo-observed accepted species percentage points per unit of mean of prior annual log(1+PCNL), original island×year fixed effects, photo observation-effort control, islands total-weight-equal, island-cluster nominal 95% CIs.

| Direct source phenotype/function documentation | beta pp | 95% CI | nominal p |
| --- | ---: | --- | ---: |
| Numerically night-higher external reproductive species (27 actual photo names) | +0.0100 | [−0.0340,+0.0541] | .655 |
| Numerically day-higher external reproductive species (25) | −0.0106 | [−0.0581,+0.0369] | .661 |
| **Paired night minus day on same source denominators** | **+0.0206** | **[−0.0438,+0.0850]** | **.529** |
| Day/night treatment hours adjusted sign comparison | −0.0244 | [−0.0782,+0.0293] | .373 |
| Source plant flowering period nighttime vs daytime | +0.01985 | [−0.0188,+0.0585] | .314 |
| Continuous published reproductive night-share centered at 0.5, same photo denom | −0.00387 | [−0.01576,+0.00803] | .524 |
| Continuous night-share normalized by day/night hours | −0.00111 | [−0.00621,+0.00399] | .670 |

### Stronger same-cell source compatibility comparison

Within the source photo cell×year with **both external night- and day-numerically-stronger plants recorded**, hold each original island and observation year fixed, require ≥2 within-island-year independent photographed cells with light range ≥0.02 and at least 25 independent supported islands.

Among the **2,373 qualified cell×years / 39 islands**, the fraction of external night-higher among photographed experimentally characterized plants has island-year/effort-adjusted historical-light coefficient:

**beta −1.17476 pp per prior mean log(1+PCNL), 95% CI [−2.47017,+0.12066], nominal p=.0742**.

This is a **selected-co-observation-only descriptive trend**, NOT statistically resolved, multiplicity-adjusted discovery, actual functional buffering, nocturnal moth decline or causal light response. Small effective independent geography (39 islands), source-study selection and unknown local species status severely limit inference.

## What was established and what was not

- Author-independent reproduction comparisons can be attached to real worldwide island flora at nontrivial scale, though the actually experimentally characterized plant species are only 52.
- **No robust global negative ALAN slope specific to plants with stronger documented nighttime reproductive response.** Broad +0.0206 and continuous −0.0039 estimates both include zero.
- Conditional same-cell −1.175 trend is a potential *study sample* pattern only, and must not be promoted by ignoring the source-wide null, hour-normalization, multiple testing and post-hoc selection.
- Flower morphology (previous PR #45) had conditional negative shape-light relation only before genus adjustment and source-quality restriction. This independent functional tier does not establish those effects either. Neither surface proves the actual environment's pollinators or plants are increasing or declining.
- Previous old GloBI −1.902 and +0.750 source-label reversal remains a different same-year VIIRS island-minus-mainland model; do not combine numerical betas from different source cohorts and exposure definitions.

## Next defensible mechanistic step

Instead of significance-hunting the conditional 39-island trend, obtain more **site-specific nighttime pollen deposition, reproductive fruit/seed outcomes**, night/day floral visitors and independent **selfing ability** with genuine wild-native/island-establishment status; repeat annual illumination and pollination assessments on the *same* source populations. Use quantitative functional redundancy (effective daytime + nighttime deposition) and independently measured predator effects as mediators only when they share true location/time support.

## Source artifacts

- PR #46 branch only; upstream PR #43–45 and main unchanged.
- Exact original author Git blobs verified, all 16 original 93m-derived island species source groups and historic PCNL signed.
- Source runner: scripts/analysis/alan_real_daynight_reproductive_PCNL_562_20261113.py.
- Synthetic source guardrails: tests/test_alan_real_daynight_reproductive_PCNL_562_20261113.py.
- GitHub workflow: .github/workflows/alan-562island-daynight-experimental-PCNL.yml.
- Public output aggregates only: no original photo/species-location, local experimental field coordinates or restricted trait cells.
