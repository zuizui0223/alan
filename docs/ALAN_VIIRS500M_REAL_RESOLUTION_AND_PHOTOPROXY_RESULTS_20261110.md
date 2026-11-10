# ALAN actual nighttime-light spatial resolution: 2024 VNL 500m, 2024 PCNL ~1km, and same-image 2x2 control

2026-11-10 | **Genuine original satellite pixels + original 2025 source photographic site IDs, SUCCESS.** No true moth behaviour/pollination or causal flower change inferred.

## Three distinct verified executions

1. [Actual fine 500m 2024 global source](https://github.com/zuizui0223/alan/actions/runs/42769904996), actual original GeoTIFF/data receipt artifact 13740914785. Real source: OpenGeoHub Zenodo 17294744 annual 2024 VIIRS VNL V2.2 15-arcsecond file, MD5 0b976892d40c2d0d99177c727289eaf2, native Int16 with GDAL SCALE=10. Publisher x10 value multiplication and GDAL SCALE=10 cancel when reading true raw pixels with rasterio. This was explicitly checked before numerical processing. Genuine original GSHHG 2.3.7 L1 landmass polygon minus lake L2.
2. [Within-original-island VIIRS500m photo-source matched contrast](https://github.com/zuizui0223/alan/actions/runs/42770395083), artifact 13741546816.
3. [Same-satellite source 500m versus 2x2 ~1km real-pixel comparison](https://github.com/zuizui0223/alan/actions/runs/42770794334), artifact 13741871816. Four authentic fine pixels are combined ONLY when all four are inside the very same 5km geodesic buffer, within the same signed island land polygon, valid and excluding lake, sea and nodata. This isolates source grid size/aggregation sensitivity without switching satellite products.

## Source PCNL versus source VIIRS, matched actual original 2025 photo sites

The fixed original 2025 photo-site pilot had 4,729 preselected sites across 1,758 source islands; 1 antimeridian original site was fail-closed before the saved source table. The source-audited valid comparison population is 4,728 source sites on 1,757 islands.

| Source exposure | Prior-year ~1km PCNL | Prior-year ~500m VIIRS VNL |
| --- | ---: | ---: |
| Genuinely pixel-supported original photographed sites | 3,543 | 3,513 |
| Supported original independent GSHHG islands | 970 | 957 |
| Median fraction of 5km island land classed dark at source-native zero | 70.09% | 65.90% |
| Median largest connected dark component / all dark land | 97.06% | 98.96% |
| Valid source sites with no zero-value dark pixels | 125 | 329 |

**3,497 same original photo sites on 951 islands are valid in both products:** Spearman source-dark fraction rho 0.9521. The within-identical-photo-site median fine minus coarse zero-dark fraction is **−1.09 percentage points**, 5th to 95th percentile −25.96 to +10.69 points. **141** same-photo sites had earlier PCNL connected dark component ≥0.95 but fine VIIRS value <0.90. These are source products with different preprocessing and photo eligibility; do not label this a pure resolution effect.

Threshold sensitivity is large: PCNL native threshold zero and one yield the same dark classification everywhere in its original valid sites, while only **22.17%** of true fine VIIRS sites have identical zero and 1 nW dark fractions. Fine VIIRS fraction dark median becomes **94.64% at 1 nW** and **99.66% at 3 nW**, making conditional dark-patch connectivity naturally saturate toward one. Thresholds are selected as sensitivity definitions, not calibrated ecological darkness at flowers.

[Committed source-verified numerical result](../results/alan_viirs2024_500m_vs_pcnl1km_darkrefugia_verified_20261110.json).

## Genuine same-image resolution-only control (NO change of VIIRS source)

A stricter all-four-land all-valid quartet selection yielded **810 independent islands and 2,756 original photo sites** with a viable 500m vs aligned ~1km comparison. At a raw VIIRS threshold of zero (both grids calculated using the exact same source pixel values):

| Median over the EXACT same 2,756 source sites | Original 500m | 2x2 1km strict dark-majority |
| --- | ---: | ---: |
| Dark land fraction | **76.79%** | **73.81%** |
| Largest connected dark fraction conditional on real dark pixels | **99.78%** | **100.00%** |

An alternative 2x2 pixel-mean light threshold yields **67.86%** median dark-land fraction — meaning the aggregation *rule* has a material role alongside pixel size. **72** authentic same-source sites have 1km strict-majority dark connectivity ≥0.95 but original 500m connectivity <0.90. The paired fine minus coarse dark fraction median is **+1.79 percentage points** at threshold zero. The paired fine minus coarse conditional connectivity median is zero; the 5th–95th quantile span is about **−0.091 to +0.125**.

This is genuine support for local grid-dependent classifications, not evidence that dark habitat is ecologically fragmented everywhere. Strict all-land quartet selection excludes coastal/small-island sites and narrows the population relative to the original 957 islands / 3513 fine-source sites. At thresholds 1 and 3 most land is classified dark, reducing the test's information content.

[Committed exact genuine same-raster resolution result](../results/alan_same2024viirs500m_vs_1km_real_quartet_resolution_verified_20261110.json).

## Biological source-bias control, not a pollination effect

Within the **2024 500m VIIRS source** and same 2025 island, after dark-fraction delta <=.10 and local log1p(mean radiance) delta <=.20, comparing sites whose conditional dark patch connectivity differs by >=.10:

| Fine 500m threshold | All original island site-pairs | Site centers >=10km apart (5km footprints do not overlap) | Higher-minus-lower photo-recorded GloBI documented moth-visitor-plant fraction (95% island-level nominal CI) |
| --- | --- | --- | --- |
| 0 nW/cm2/sr | 212 islands / 479 pairs | 100 islands / 188 pairs | **+1.03 pp** [−3.26,+5.31] |
| 1 nW/cm2/sr | 120 islands / 298 pairs | 64 islands / 126 pairs | **−5.07 pp** [−10.42,+0.27] |
| 3 nW/cm2/sr | 94 islands / 250 pairs | 55 islands / 118 pairs | **−4.12 pp** [−9.85,+1.61] |

No consistent source-plant positive relationship survives these threshold and non-overlap checks. All confidence intervals include zero; no correction for the multiple exploratory specifications has been applied. Plant observation effort/recorded richness controls remain nondefinitive. GloBI recorded moth-visited plants are not local observed moths, obligate or effective nocturnal pollination, insect population abundance, or wild native plant persistence.

[Committed observed source photo test](../results/alan_viirs2024_500m_source_photo_proxy_matched_verified_20261110.json).

## Decision

**Keep 500m VIIRS VNL V2.2 as the reproducible global island spatial source, PCNL as historical coarse 1992–2024 exposure, and the same-source 2x2 audit as the resolution-uncertainty test.** The fact that different thresholds change spatial classification and photographed source-proxy conclusions warns against declaring a robust ecological effect. No sub-100m streetlamp attraction can be inferred from these rasters. SDGSAT-1 provides 10m PAN/40m RGB night imaging but requires local scene access/approval; not a global substitute. Further biological testing requires independent same-site nocturnal insects, visits, pollen deposition or seed/fruit measurements.
