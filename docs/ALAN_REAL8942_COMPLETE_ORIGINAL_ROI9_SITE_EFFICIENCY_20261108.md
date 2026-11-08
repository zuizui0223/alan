# ALAN — genuinely new original RGB observations: complete 8,942-photo ROI9 source measurement (2026-11-08)

**Completed empirical data acquisition and independent source-only re-verification.**

[Original 7,662-photo 15-worker run](https://github.com/zuizui0223/alan/actions/runs/42456363487) is marked **FAILURE** because of a **shell error in the log-printing line, AFTER each original 494–512-photo FCP inference unit finished**. All fifteen workers nevertheless uploaded their authentic original measurement terminal photo CSVs plus all four per-128-photo original SHA receipts. Crucially, **no new download or re-inference was used** to recover these results.

[ALAN independent source and actual-photo receipt-recovery run #42465405115](https://github.com/zuizui0223/alan/actions/runs/42465405115) is **SUCCESS**. It rehydrated all 15 saved worker datasets, recreated and byte-compared original source manifests, verified 60 source shard SHA/terminal CSV pairs and their unchanged location-blind FCP YOLO11/EfficientSAM model blob, then combined with previously verified original source batches 1–5. **All 8,942 photo IDs, 8,942 observation UUIDs, original island IDs and photo-cell×year source IDs matched the precolour, signed 35-batch source roster**.

## What actually grew

| Real source photo result | Actual |
| --- | ---: |
| **New independent original observation UUIDs / original image attempts** | **8,942** |
| **Original numeric nine-component floral-region ROI9 successes** | **3,597 (40.2%)** |
| Technical unknowns including original-image retrieval failures | **5,345** |
| Permanent original-image acquisition failures | **1** |
| Distinct source GSHHG islands | **369** |
| Distinct original photo-cell×observation-year units | **4,776** |
| Site×years with at least one actual successful numeric floral ROI9 | **2,517** |
| Site×years with source opportunity ≥5 original photos + ≥3 observed WCVP accepted taxa | **303** |
| **Site×years actually passing ≥5 photo attempts, ≥3 successful numeric ROI9 photos, and ≥3 accepted taxa with successful original ROI9** | **156** |
| Distinct GSHHG islands in the actual 156 strict site×years | **12** |

The 156 strictly eligible *new acquisition-arm* source photo cell-years break down by original observation year: **2019: 4; 2021: 15; 2023: 34; 2024: 44; 2025: 59**. Since they are concentrated in only 12 islands, **they do not by themselves solve global island representativeness**, and must not be silently added as randomized independent support to the original 710 island photo-cell×years / five-year global historical-main regression.

The measured response remains the **numeric nine-palette mass of a photographed flower-region ROI**, not verified biologically petal-only colour nor a census of plant species community abundance. Technical unknowns were never assigned a white or coloured label.

## Evidence of acquisition design efficiency, now verified on the same REAL photos

The separately source-preselected 512 observation-UUID **site-year densification** records are a subset of the now-measured 8,942 new independent original source photos, with original source SHA signed before their colour outcomes were used. The original first 1,280 photo IDs were held unchanged.

| Strict photographic site comparison | First 1,280 + next 512 geographically broad | First 1,280 + 512 targeted within real site-year cells |
| --- | ---: | ---: |
| Distinct original photo attempts | 1,792 | 1,792 |
| Original numeric ROI9 successes | 675 | 691 |
| Precolour source 5-photo/3-observed-taxon opportunity sites | 6 | 60 |
| **Actual 5-photo/3-successful-ROI/3-successful-taxon sites** | **1** | **33** |
| Distinct sampled islands in combined first 1,792 | 369 | 369 |

The source-targeted comparison yields **33 times as many usable source photograph cell×years** with essentially equal original photo acquisition workload, and only 16 more measured numeric ROI9 photos. The improvement arises largely from **clustering existing observation opportunities into comparable units**, not from a miraculous image inference success-rate change. The designs are purposive and overlap in their first 1,280 source photos; they are **not independent randomized trials**. Conditional environmental/geographic sampling and original-image detection still require modeling before drawing flora ecology conclusions.

## Actual model, failure mode, correction and verification

1. The **original full 93,209,351-photo iNaturalist index** had already been partitioned into immutable source metadata tiers. The tier06 record supplied a precolour schedule of **8,942 actual new observations** with same-year existing VIIRS/WorldPop support, one source photo ID per observation.
2. Original source-specific S3 photos were acquired and the unchanged YOLO11/EfficientSAM FCP original ROI-v4 evaluator run in bounded warm model workers, 1–5 initial groups and 15 additional paired groups. Original photo bytes were not persisted; each worker stored original photo-level numeric and unknown outcomes and SHA-signed source receipts. The 15-worker run log was incorrectly terminated by a shell variable in a final echo, **after photo FCP inference**. It therefore must **not** be called a successful Actions run, even though all actual original photo results were present.
3. Recovery used **only completed original photo CSVs**, not downloaded RGB pixels. It recreated every original source shard from the signed source roster, verified per-shard photo ID/observation UUID and photo-result SHA, and required the identical frozen FCP source Git blob `b1d7091db44682559bbe4b48a9f3b8db648501a6`.
4. The successful data-only run [#42465405115](https://github.com/zuizui0223/alan/actions/runs/42465405115) validated every signed source image outcome, all 8,942 unique original photo IDs and 8,942 distinct observation UUIDs, and the **numerical nine-palette success/unknown decomposition 3,597 + 5,345 = 8,942**. Per source year, photo cell and originally accepted WCVP plant taxa, it counted only **actual** 5-attempt/3-numeric-ROI/3-accepted-taxon strict site×years.
5. The original photos' acquisition arm is explicitly island-prioritized. The older historical-main island×ALAN regression results were **not refitted** and MOIS remains for moisture×flower colour.

## Reproducibility and source retention

- [Original full 8,942 precolour photo source SHA schedule](https://github.com/zuizui0223/alan/actions/runs/42447983019): `e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd`.
- [Successful data-only full 8,942-source ROI9 audit run](https://github.com/zuizui0223/alan/actions/runs/42465405115), original result artifact **`alan-8942-real-original-island-ROI9-full-source-dataonly-recovery`** ID **13611338201**; final source receipt JSON SHA256 **`0043a4b2abedc613136d6a6544b3b3a70282fb707d5fcc760e0cc5a1111e06d2`**.
- Full signed source receipt and a per-35-batch tabulation in the artifact; **[compact GitHub result JSON](../results/alan_real8942_original_nineROI_verified_summary_20261108.json)**.
- Reproduction scripts: `scripts/acquisition/alan_bulk_7662_original_roi9_source_20261108.py`, `scripts/analysis/alan_collect_full8942_island_original_roi9_20261108.py`, QA tests and the `.github/workflows/alan-8942-originalROI9-dataonly-recovery-audit.yml` workflow.
- Original photo and prior MOIS signed artifacts have finite GitHub Actions retention. No permanent new imagery archive was created and no rights to redistribute image bytes are implied.

**Next scientific objective:** choose additional source photo opportunities on islands that have enough accepted taxa and photo observations for independent strictly measured *site-years*, while retaining a geographically spread reference cohort. The newly source-verified tier07 27,507 original island observation candidates form the next acquisition population; they have not yet been RGB-classified in this source study.
