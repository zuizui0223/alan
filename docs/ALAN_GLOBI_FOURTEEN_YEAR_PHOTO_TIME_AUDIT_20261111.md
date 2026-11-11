# GloBI × artificial night lights: long 2012–2025 botanical-photo time audit (2026-11-11)

**Executed real source audit. The new 14-year estimate is a photo-observed plant-species composition trend, NOT a measured 14-year ALAN interaction.** Existing `fd` and `island` Chapter 1 were not changed.

## Research extension

The previous confirmed empirical ALAN series used 2019–2025 source photos and genuinely same-year (but sampled only Jan/Apr/Jul) VIIRS for 7 calendar years. The user asked for a longer time series. We returned to the original GSHHG/WCVP-signed image observation archive and, using ALL 61,967 formerly acquired positive GloBI Lepidoptera–plant claims, reconstructed the **2012–2025 photo-species composition** without requiring an `fd` morphological phenotype.

- Full original historical_main photos for 2012–2025: **649,694**, including **63,715 from 2012–2018**. Source SHA256 `f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c`.
- GSHHG 2.3.7 L1/ L2 point-in-land geometry, with inland water exclusion and correct physical sibling landmass IDs, recomputed for early years and all of 2020. Private exact WCVP accepted species and photo-cohort source retained. No original image downloads.
- **332,576 observed photo-site×year rows** across the 14-year combined source; **54,055** have at least 3 independent photo opportunities and 2 accepted photographed plant species.
- **10,478** broad-source photo sites have at least two observation years, including **1,316 island sites**.
- Direct same-calendar-year and same-5° region island/mainland support: **14,590 photo-site×years (1,338 island)**, 85 distinct 5° regions.

### Empirical observational photo-composition slope

Endpoint = percent of actually photographed accepted plant species with a known positive GloBI moth-family flower interaction *somewhere* in its source archive. Not the count of local moth visits, actual pollen transfer, native plant richness, absence or population decline. Equalize total island/mainland weights within every observed year×geographic block. Control photo opportunity, observed accepted species richness and (in sensitivity) photo month and hemisphere; geographic 5° cluster-robust exploratory uncertainty.

| Original photo scope | Island-minus-mainland yearly change in photographed GloBI-moth-positive plant fraction | Nominal 95% CI | p |
| --- | ---: | --- | ---: |
| **2012–2025 full fourteen calendar years** | **−0.0068 pp/year** | **[−0.735,+0.721]** | **.985** |
| 14 years, photo-season and hemisphere adjustment | −0.0459 pp/year | [−0.838,+0.746] | .909 |
| 14 years, true photo-site fixed effects (7,635 site×years; 2,479 sites) | −0.516 pp/year | [−1.121,+0.089] | .092 |
| 2013–2024 photo scope (annual VIIRS **NOT actually joined**) | −0.330 pp/year | [−1.207,+0.547] | .455 |
| 2017–2025 better-supported photo geography | −0.0178 pp/year | [−0.811,+0.775] | .965 |

A **post-hoc** 2017–2025 same-site FE sensitivity gives −0.870 pp/year [−1.628,−0.112], nominal p=.0256, but **loses support when iNaturalist-derived GloBI interactions are removed**: −0.403 [−1.010,+0.204], p=.187. Do not interpret this selected result as a verified long-term loss of moth-pollinated island plants.

### Independent early vs late sites

Among **1,848 photo sites** represented in both 2012–2018 and 2019–2025 (island 289/mainland 1,559), only **40 island and 537 mainland sites in 16 geographical 5° blocks** have independent same-region two-class temporal common support. Compare each site's early vs late observed positive-moth photo-species fraction, then its within-region mean delta difference. Equal-block mean **−4.324 pp**, block-bootstrap 95% CI **[−17.092,+8.110]**, exact whole-block sign-swap p **.531**. No robust relative loss.

### Pre-2017 source weakness is substantial

Under the 3-photo/2-accepted species gate, island source points in the *same-year same-geographic block as mainland points*: 2012 **1**, 2013 **7**, 2014 **0**, 2015 **4**, 2016 **4**, 2017 **13**, 2018 **40**, rising substantially from 2019 onwards. No imputed absent island-year controls; 2014 contributes no directly comparable island data. Unbalanced observation intensity rises from **3,222 historical_main original photos in 2012 to 125,033 in 2025**.

### 2020 exact source reconciliation

Old 2020 VIIRS partial geographic source contained **25,327 sites**; recomputing the actual canonical 2020 image site population gives **31,102** source site-years. In the **25,290** sites with exact target ID overlap, 233 had changed photo/opportunity or species/interaction counts after exact source re-construction; the original and rebuilt populations must not be claimed byte-identical. The 14-year basic photo-composition year slope is −0.0068 on full reassembled 2020, −0.134 on the former partial 2020, and +0.033 if 233 differing sites are dropped. None establishes a biological long-term trend.

## What has NOT been estimated

An ecological **2012–2024 annual-nightlight × island × GloBI-positive plant** interaction has NOT been executed: source-matched official annual VIIRS radiances, geospatial masks, `cf_cvg` cloud-free QA and source release provenance for original 2012–2018 photo sites remain unavailable in this source run. Likewise, no 1992–2011 DMSP-OLS source light exposure was attached. 1992–2011 historical_main original photos number **17,126** only, so the old era is also sharply under-supported botanically.

The appropriate official EOG annual VIIRS sources are [VNL V2.1](https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_ANNUAL_V21) for **2013–2021** and [VNL V2.2](https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_ANNUAL_V22) for **2022–2024**; 2012 has partial/unlike processing in the original catalog. Do NOT treat 2025 as the same annual source. A strict input checker and `island × log(1+annual radiance) × calendar year` model are included in the local reproducibility ZIP; it **fails closed** on missing actual radiance, source-year mismatch, wrong version or wrong bands. These tests passed on synthetic input, but no annual source/light effect was generated.

**Scientific interpretation:** The old 2019–2025 source-light spatial association is not an established biological temporal trend. The newly extended 14-year *photographic* source does not yield a consistent island-specific monotonic trend in the documented moth-positive plant fraction. It cannot by itself establish an ALAN effect, moth decline, naturalization failure, gecko mechanism or species replacement.

## Reproducibility and rights

Public aggregate-only outputs:
- `results/alan_GloBI_long_photo_2012_2025_source_verified_20261111.json`.
- User-facing ZIP `ALAN_GloBI_14year_photo_temporal_20261111.zip`: exact source-geography code for the earlier years, actual concatenated photo-year model, 2020 source reconciliation, secure annual source gates, 9 passed synthetic tests and aggregate-only results. **The raw species×site/year/photo coordinates are not in either output.**

The original 2019–2025 same-year VIIRS coefficients, prior source hashes and original GloBI cache remain frozen; this is an independent extra long-time-series observational arm. Draft PR must remain unmerged until official annual radiance and robust source-independent visitor evidence are established.
