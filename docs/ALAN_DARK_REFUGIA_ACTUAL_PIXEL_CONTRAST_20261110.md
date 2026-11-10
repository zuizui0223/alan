# ALAN — actual dark-refuge spatial measurements and independent conditional support, 2026-11-10

**Status: source-verified original satellite measurements SUCCESS. NOT a biological ALAN–moth or ALAN–plant effect.**

Two successful GitHub Actions:
- [Original 2024 PCNL pixels + original GSHHG landmask measurement](https://github.com/zuizui0223/alan/actions/runs/42764291562) — artifact 13738439540; real gzipped 2025 site-level 5-km satellite-pixel exposure table, original table SHA256 32327abb37ee85f3d1d9e7446d9d742fb5adf1457dab60fae1aae791c29f7837.
- [Observed source support, frozen thresholds and matching diagnostics](https://github.com/zuizui0223/alan/actions/runs/42765060768) — artifact 13739591784. Numeric JSON committed at ../results/alan_2025_dark_refugia_conditional_spatial_support_20261110_verified.json.

## Genuine source support

2025 observed plant-photo same-island cell-years: 4,729 SHA-sorted pilot sites on **1,758 independent source GSHHG islands**, at most eight actual photo sites per island, chosen without reading any biological response. One antimeridian site fails closed (unsupported wrapped raster). 4,728 yield satellite neighbourhood diagnostics, of which **3,543 independent 0.05-degree photo sites on 970 islands** meet original 2024 pre-observation PCNL five-km landmask coverage gates (>=16 source land pixels, >=80% valid coverage).

For native PCNL dark threshold = 0 (not lux, not verified no ground lighting):
- Median measured dark-land fraction: **0.700862** of land in the 5 km area; 5th–95th percentile 0.03036–1.0000.
- Among photo sites with at least one dark pixel: median fraction of dark land within the largest four-connected dark patch: **0.970584**; 5th percentile 0.33338 and 25th percentile 0.68424.
- 125 valid-photo sites have **no** dark PCNL pixels within that 5 km landmask. For them the fraction within a largest dark component is undefined; nearest-dark distance is censored at radius, not zero.
- At threshold PCNL=1, **all 3,543 sites have identical dark-land fractions** to threshold=0. These thresholds are redundant for this original, coarse, quantized product, not independent robustness checks.
- At threshold PCNL=3, dark-land fraction median **0.776095** and conditional dark patch continuity median **0.979182**.

## The critical independent-variation test

Of 970 source-qualified islands, **659** have >=2 distinct source site footprints for within-island comparison. At threshold PCNL=0 there are 7,627 pairs of actual dark-containing photo sites, enumerated within same island and same original observation year only, without treating threshold repeats as replicas.

Prespecified exploratory environmental common-support criteria: difference in *dark-land fraction* <=0.10, difference in local log1p(mean PCNL light) <=0.20, but difference in conditional largest-dark-patch connectivity >=0.10 (all measured, not fitted). Result: **544 real within-island site pairs across 235 source islands**.

Without loglight matching there are 1,039 pairs / 308 islands meeting dark fraction and connectivity comparison; requiring a stronger >=0.20 difference in connectivity gives 666 pairs / 224 islands without additional light matching. At threshold 3, the full matched-light criterion gives 549 pairs / 230 islands.

Thus the needed landscape **spatial contrast is actually observable** on a moderate number of original photographed islands, in spite of the high median connectivity. It does **not** prove an ecological effect, because no true moth-availability, flower-visit, pollen deposition or fruiting response has been linked at these matched sites. The plant-photo observation site is not a surveyed floristic quadrat.

## Next questions and limits

- First test whether a truly independent moth behavioural outcome is available at exactly these sites, or whether site-level detection effort/source selection makes any biological identification invalid. Do NOT use GloBI source-positive moth-visited plant ratio as effective moth abundance or a true negative.
- Test 5- and 10-km landmask source resolution/scale robustness, holding the same island cohort and fixed prior-year exposure. Quantify spatial overlap of neighbourhood footprints within island rather than counting highly overlapping sites as independent ecological replicates.
- The four-connected PCNL pixels use a coarse approximate ~1 km source pixel; these **cannot identify streetlight attraction within 0–100 metres**. They measure *landscape-scale* dark mosaic only, not individual streetlight attraction, moth flight corridors or habitat suitability.
- Fragmentation from the same thresholded light raster is necessarily correlated with dark fraction. Matching on the latter, mean light and source geography is necessary but not sufficient for causal claims; land cover, urbanization, elevation, botanical photo effort and moth-trap design would remain major confounders.
- All inferences are post-exploration and source-chosen, not confirmatory. No image re-inference or earlier flower-white model changes have been made.

**Conclusion:** the geospatial **predictor** is constructed and its independent within-island contrast is measurable (235 islands). The ecological **outcome** and mechanism remain unverified; no nightlight-induced positive or negative moth effect has been estimated.
