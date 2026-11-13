# ALAN: 562-island true-source moth-like flower morphology × historical nightlight

**2026-11-13; real signed-data rerun complete, not a biological pollination experiment.**

Reproducibility: [successful GitHub Actions #43340299255](https://github.com/zuizui0223/alan/actions/runs/43340299255), [aggregate-only JSON artifact #13999495995](https://github.com/zuizui0223/alan/actions/runs/43340299255/artifacts/13999495995), [Draft PR #45](https://github.com/zuizui0223/alan/pull/45).

## Exactly matched source and estimand

- Original 93,209,351 photo-ID universe; 16 original accepted WCVP species × island/site/year SQLite source groups. Fixed actual macro sample: **562 islands, 77,647 cell×observation years, 2,012,925 species×site×years, 14,330 distinct species**. Actual observation years 2019, 2021, 2023, 2024 and 2025. These are not 93 million flowering-specimen images.
- Predictor: authenticated Li–Cao corrected annual PCNL, 1992 to **year preceding each plant photograph**; average annual log1p values used as long-term previous light exposure. Not ground-level lux or a long-term population trend.
- Species-DIRECT high/medium shape+color+depth from independently SHA-verified Island 2026-09-08 134,268-row trait archive. No model/genus proxy. The original high-only quality restriction is re-fitted separately. Source rights protect species–trait rows; outputs aggregate only.
- White + tubular / salverform / funnel-trumpet is a **moth-LIKE morphology**, not demonstrated local moth pollination; strict subset additionally needs direct species evidence of deep flower tube. Missing original trait means UNKNOWN, not noncandidate.
- Original contemporaneous species GloBI prior-PCNL source coefficient **+0.4572282788** independently reproduced on the exact same original site data as a provenance gate. Do not compare that directly with the distinct island-minus-mainland same-year VIIRS old −1.902 coefficient.

## Species support and documentation groups

| Direct botanical trait quality | Both color+form known | Broad white+tubular-like | Strict white+tubular+deep |
|---|---:|---:|---:|
| High+medium | **1,691 species** | **257 species** | **13 species** |
| High-only | **811 species** | **114 species** | **3 species** |

Broad high+medium candidates were actually photographed in **444 islands / 21,984 cell-years**, high-only **348 islands / 12,804 cell-years**. Strict high+medium documented plants occurred across 88 islands / 1,145 photo-cell years, but their effective nocturnal pollination is unobserved.

Within the new 14,330-species frame (NOT the earlier 8,581-species population): old-only GloBI 784 taxa with 150 paired direct trait records and **8 broad candidates**; old+new recorded 637 with 109 known and **9** broad candidates; new-only 163 with 40 known/**5** broad candidates; neither source recorded 12,746 with 1,392 known/**235** broad candidates. Exploratory two-sided Fisher for broad morphology among known records, old-only 8/150 vs retained 9/109: **odds ratio 0.626, p=.447**. No old-only candidate enrichment. High-only: 4/99 vs 3/66, p=1.0.

## Verified empirical historic-light coefficients

All β are photo-recorded species composition percentage points per original mean prior log1p PCNL, not species disappearance or any actual pollen transfer. Source main model: island×observed-year fixed effects; log original observation effort; each independent island has equal total weight; island-clustered nominal 95% interval.

| Source-photo botanical endpoint | High+medium β [95%CI] (p) | High-only β [95%CI] (p) |
|---|---|---|
| **White+tubular-like** / all photographed species | **−0.0834 [−0.1693,+0.0024] (.0567)** | **−0.0106 [−0.0689,+0.0478] (.722)** |
| White flower alone | −0.2823 [−0.5231,−0.0415] (.0217) | −0.0090 [−0.1972,+0.1793] (.925) |
| Tubular-like form alone | +0.3831 [+0.1284,+0.6378] (.0033) | +0.1433 [−0.0578,+0.3445] (.162) |
| **Florally recorded FORM at all** | **+0.5419 [+0.1497,+0.9340] (.0069)** | +0.3516 [−0.0249,+0.7280] (.067) |
| Both shape+color documented | −0.0997 [−0.3553,+0.1560] (.444) | −0.1130 [−0.3418,+0.1158] (.333) |
| White+tubular+deep | +0.0052 [−0.0121,+0.0225] (.553) | +0.0049 [−0.0043,+0.0141] (.293) |

**Trait ascertainment warning:** the shape-recorded indicator itself is positively associated with nightlight (nominal p=.0069 in high+medium). The apparent white-only negative and tubular-only positive slopes fail independent high-quality trait restriction. Post-hoc many correlated outcomes are not independent discoveries.

### Primary stricter comparison: only fully measured color+form species

After ≥2 independent morphology-known recorded species per cell-year and genuine contemporaneous island-year light support, high+medium includes **51,524 photo cells across 433 islands**; broad candidates occur at 19,627 cells in 351 of these islands. The outcome is candidate versus measured noncandidate, not an undocumented species called negative.

| Model | High+medium β [95% CI] (p) | High-only β [95% CI] (p) |
|---|---|---|
| Known-only, no genus adjustment | **−0.7905 [−1.2724,−0.3086] (.00136)** | −0.5113 [−1.0680,+0.0453] (.0717) |
| **Known-only, genus fixed effects** | **−0.0228 [−0.1074,+0.0618] (.596)** | **−0.00388 [−0.0514,+0.0437] (.873)** |
| Only genera containing both morph classes, genus FE | −0.3663 [−0.9627,+0.2301] (.227), 53 genera | −0.7943 [−1.7082,+0.1196] (.087), 12 genera |

The apparent negative association in original studied species composition is **greatly attenuated once plant genus is controlled**. Do not call the model difference a tested biological mediation effect: genus can absorb biogeography, documented species pool, observer preferences, evolutionary lineage and shared environment simultaneously. The higher-quality source independently weakens the raw association.

Old-only GloBI-positive broad morphology is **eight species** in this 562-island frame, with old-only morphology source photo gradient +0.00684 [−0.01745,+0.03113], p=.581. It does **not** recreate the earlier GloBI documentation label source-negative contrast. The two exposure regimes, source cohort and response denominators differ and their β magnitudes cannot be subtracted.

## Scientific claim boundary and next step

**Current conclusion:** There is no robust supported island light association *specific to morphologically moth-like plants* after genus and trait-quality controls. This does NOT prove ALAN harmless to effective nocturnal pollination. The observed source is neither a wild/native island flora census, real at-flower pollinator census, true effective nocturnal pollen-deposition assay nor a repeated-plant-population trend. Most candidate evidence remains partially measured; deep floral tubes are especially under-documented.

Further sign hunting across area/latitude/assumed absence is inappropriate. Next mechanistic test should obtain independent genuine nocturnal pollen-transfer/day-night exclusion evidence and island-wild status for a predeclared set of species, then compare net plant reproductive outcomes and actual lights/dark landscape arrangement; retaining true unknown traits as missing.

Source scripts: scripts/analysis/alan_real93m_direct_flower_mothlike_PCNL_20261113.py; source-provenance helper scripts/analysis/alan_real_genus_prevalence_bias_20261110.py; tests/test_alan_real93m_direct_flower_mothlike_PCNL_20261113.py; .github/workflows/alan-562island-93m-direct-floral-morph-PCNL.yml. Original source trait identities are NOT published in these outputs. Earlier Draft PR #43/#44 and main remain unchanged.
