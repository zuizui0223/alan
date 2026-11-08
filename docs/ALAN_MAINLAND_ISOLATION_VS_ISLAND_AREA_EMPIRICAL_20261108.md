# ALAN 2026-11-08 — Does mainland isolation explain the island-area × nightlight gradient?

**Status:** Real original-photo five-year analysis executed in ALAN, all job steps **SUCCESS**. [Execution #42440394223](https://github.com/zuizui0223/alan/actions/runs/42440394223) · [Inspect the signed aggregate JSON summary](../results/alan_isolation_area_light_20261108_source_verified_summary.json). This is an exploratory, observational, not causal, result.

## Research question and precise measurement

The previous ALAN study found that **larger islands have a less negative observed association of nightlight with photographic floral-region whiteness**. Because small island area and geographic isolation are different properties, the new question is:

> Is the **interaction of island area and artificial nightlight** still observed when the interaction of **distance to continental coastline and nightlight** is included in the *same* model?

Using the same source SHA-gated GSHHG v2.3.7 L1 coastline from the previous 5-year studies, calculate **minimum WGS84 geodesic among GSHHG original outer-coastline vertices** of each selected island and the original study's GSHHG landmasses with **area greater than 7 million km²**. Antimeridian-safe ECEF nearest-neighbour search screens three nearest mainland vertices for each original island coastline vertex. The resulting number is a **coastal-vertex approximation and an upper-bound to exact polygon-edge minimum**. It is *not* distance to nearest island, centroid-to-mainland distance, distance along a ferry route, or independently measured pollinator isolation. No new environmental or photographic values were invented.

Define `A = log10(original GSHHG island area in km²)`, `D = log10(1 + vertex-distance km)`, `L = log(1 + annual VIIRS)`. The model among previously eligible multi-site islands is:

```
white_photoROI_pp ~ L + L:A_centered + L:D_centered
                   + log_pop + BIO5 + log_rain + elevation + log_photo_attempts
                   + island fixed effects + observed-year fixed effects
```

The individual interaction slopes are **pp of photographed floral ROI white mass per unit L, per tenfold area or tenfold (1+distance)**. All preexisting 2019, 2021, 2023, 2024 and 2025 source joins and island/photo identity requirements were unchanged, preserving the strict 44-island / 650 photo cell-year and broader 95-island / 1,884 photo cell-year subsamples. Mainland records are not included in this particular within-island model. Errors are clustered by original GSHHG island, years are not separate island replicates. The analyses are both cell-year-equal and **total-weight-one-per-island**, the latter preventing highly photographed large islands from dominating.

## Observed **joint** area and coastline isolation interactions

| Cohort and observational weighting | Area × light β (95% CI) | Mainland distance × light β (95% CI) |
| --- | --- | --- |
| **Strict, island-equal (44 islands)** | **+5.384 [2.678, 8.089]** | −1.365 [−5.106, 2.376] |
| Broad, island-equal (95 islands) | **+2.939 [0.833, 5.045]** | +1.841 [−0.782, 4.463] |
| Strict, cell-year-equal | +2.217 [−0.004, 4.437] | −1.257 [−3.096, 0.583] |
| Broad, cell-year-equal | +1.016 [−0.384, 2.416] | −0.148 [−1.349, 1.053] |

Nominal exploratory cluster p for *area* in joint island-equal models: strict **0.000236**, broad **0.00674**. The isolation joint estimates' nominal p: strict **0.466**, broad **0.167**. No joint isolation interval excludes zero. The **same positive sign for area** appears across all four joint specifications, but the **cell-year-equal** confidence intervals include zero; their point estimates should not be called independently confirmed positive effects.

### Source support, collinearity and reconciliation

- Strict source **44** islands / **650** authentic cell-years, median nearest mainland coastline vertex distance **80.70 km**; broad **95** islands / **1,884** authentic cell-years, median **125.96 km**.
- Equal-island correlation of `log10 area` and `log10(1+distance)`: strict **+0.354**, broad **+0.200**. Thus the two source metrics are **not interchangeable**, though neither correlation alone proves full statistical identifiability or ecological mechanism.
- Existing **area-only** island-equal estimate is independently reproduced exactly within tolerance: strict **+5.152 [2.559, 7.744]**, broad **+3.385 [1.319, 5.450]**. **Adding continental isolation does not remove the positive area interaction** under island-equal weighting.
- The **isolation-only** model is unstable across source cohorts. For broad island-equal only, it yields nominal **+2.794 [0.035, 5.553]**, whereas the joint isolation interval includes zero. Do **not** claim an established directional isolation effect on the basis of this lone exploratory model.
- Identical original-photo/ROI/accepted taxon sources were reused without re-imputation. No newly processed photos or ROI inferences. The output is aggregate only; individual island locations or licensed original photographic records are not redistributed.

## Ecological interpretation and manuscript boundary

**Positive observed pattern retained:** Within islands, the correspondence of nightlight with photographed floral whiteness changes with island area. Even after introducing an independent high-resolution continental-distance proxy, this positive interaction persists in equal-island models on the strict and broader observed-source cohorts.

**Not established:** Geographic isolation as an independent modifier; true plant-community turnover or native white-flowered species disappearance; ALAN exposure to petals at ground level; moth abundance, effective nocturnal pollination or light-caused natural selection. A correlation with island size can also arise from land use, habitat heterogeneity, harbor/maritime development, photographed taxonomic composition, or biased photographer recruitment.

The GSHHG **landmass definition** treats only polygons larger than 7 million km² as 'mainland'. Consequently the distance estimate is **not distance to the nearest source population or nearby island**. The nearest *vertex*, not continuous coastline segments, provides an approximate upper bound. These are explicit biological/methodological limitations, not reasons to discard the measured positive area pattern.

**Next substantive test after this stage:** connect the area-related photographic colour gradient to independent source traits or actual pollinator community data on the *same sampled islands*, without treating global moth visit records as proof of local effective nocturnal pollination.

## Source verification and files

- Source-verified run: [42440394223](https://github.com/zuizui0223/alan/actions/runs/42440394223), **SUCCESS**. Artifact `alan-5year-GSHHG-island-isolation-area-nightlight`, ID **13601153247**; internal `island_isolation_area_light_aggregate_receipt.json` SHA256 `668e1ca74b84ab2538b5b4bd8d1b5c7144b2828f23e9af1b5bcb18617958a441`.
- Source GSHHG archived ZIP SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`. All other source SHA fields are retained in the archived full JSON.
- Reproducible script: `scripts/analysis/alan_isolation_area_light_20261108.py`; geographic/fit synthetic unit tests: `tests/test_alan_isolation_area_light_20261108.py`; standalone source-verified Actions job: `.github/workflows/alan-isolation-area-nightlight.yml`.
- The historic MOIS raw artifact dependency and its eventual expiry have **not** been resolved by this statistical analysis. The old results and photographic data are left in place; MOIS remains moisture-only for future research.
