# ALAN — tier07 first 1,024 independent island observations: recovered original ROI9 empirical result (2026-11-09)

**Actual original image inference and source verification complete.** The initial [GitHub Actions #42466504871](https://github.com/zuizui0223/alan/actions/runs/42466504871) stopped at a 70-minute **timeout** while the second 512-image subgroup was being processed. Its uploaded source-signed artifact contained the original photographic ROI9 terminal records and SHA receipts for **7 completed 128-photo work units (896 genuine independently observed photos)**.

The later [data-preserving recovery #42552642250](https://github.com/zuizui0223/alan/actions/runs/42552642250) completed **SUCCESS**. It regenerated the eight source photo rosters from the SHA-pinned [27,507 observation precolour source](https://github.com/zuizui0223/alan/actions/runs/42466279952), byte-compared the signed source manifest, verified the earlier seven FCP photo results with their original SHA source rosters, then **retrieved and inferred only the missing final 128 original photos**. There was *no re-download or re-inference* of the preceding 896 photos. All 8 original source ROI9 terminal receipts were checked with the exact previously frozen YOLO11/EfficientSAM FCP-v4 method.

## Actually completed original photo colour measurements

| Outcome in original source, not botanical flora prevalence | Verified result |
| --- | ---: |
| Independent *new* original observation UUIDs and original photographs | **1,024** |
| Distinct GSHHG islands with original photo source observations | **284** |
| Distinct original 0.05° photo cell × observation-year locations | **324** |
| **Actual numeric original nine-component flower-region ROI9 photographs** | **403 (39.4%)** |
| Actual technical ROI9 unknowns (not coloured or white by imputation) | **621** |
| Original-image retrieval failures | **0** |
| Original photo site×years with ≥5 source photo attempts and ≥3 accepted taxa (before pixel inference) | **76** |
| **Actual strict site×years with ≥5 attempts, ≥3 authentic numerical ROI9 successes from ≥3 accepted taxa** | **57** |
| Distinct GSHHG islands represented by these 57 actually strict sites | **29** |

The 57 real strict photo cell-years occurred in **2019: 3; 2021: 5; 2023: 13; 2024: 13; 2025: 23**. They constitute a **purposive new acquisition arm**, not an independent, probability-representative extension of the previous historical-main 710 strict island photographic cell-years.

The original source split was determined **before opening pixels**, by photo/observation IDs and source site-year opportunities, not by measured white colour:

| First 1,024 original photo acquisition part | Original photos | Authentic numeric ROI9 photos | Distinct source islands |
| --- | ---: | ---: | ---: |
| Geographical reference coverage | 256 | 112 | 256 |
| Accepted-taxon-rich within-site×year densification | 768 | 291 | 36 |
| **Total** | **1,024** | **403** | **284 distinct** |

Islands overlap between the two parts, so the separate island counts must not be added. The source-design experiment is observational and *not* proof that ALAN causes loss of white-flowered plants, that moths disappeared, or that flower colour evolved. Original source photos are not an independently confirmed census of wild native plants in anthesis. The 57 eligible cells are **genuine data-quality qualification**, not an effect-size estimate.

## Precisely verified files and SHA receipts

- [Complete 1,024-photo real FCP recovery Actions](https://github.com/zuizui0223/alan/actions/runs/42552642250): SUCCESS; artifact `alan-tier07-original1024-measured-only-last-shard-recovered-ROI9` ID **13645374703**, containing original measured 9-palette CSV, eight original 128-photo source terminal receipts, and source manifests.
- Original final numeric ROI file SHA256 `2988a8631cfb25bea57690a2a91db1e0a7f00074b8ca3ce7557aba8c09fda113`.
- [Precolour original 27,507 new observation master](https://github.com/zuizui0223/alan/actions/runs/42466279952), SHA256 `7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62`.
- Unchanged FCP evaluator original source Git blob `b1d7091db44682559bbe4b48a9f3b8db648501a6`.
- [Machine-readable, source-backed numerical result](../results/alan_tier07_first1024_actual_originalROI9_20261109_verified.json).
- Original source and photo pipeline `scripts/acquisition/alan_tier07_27507_photo_source_schedule_20261108.py`, `scripts/acquisition/alan_tier07_originalRGB_roi9_1024pilot_20261108.py`, and recovery job `.github/workflows/alan-tier07-first1024-recover-only-missing-shard007.yml`.

The remaining **26,483** original photo/observation UUIDs in the same disjoint tier07 source schedule are processing in [52 bounded warm RGB/ROI9 jobs](https://github.com/zuizui0223/alan/actions/runs/42553223750). Those remaining 26,483 image outcomes have **not** been included in this 1,024-photo verified report. The corresponding [Draft PR #13](https://github.com/zuizui0223/alan/pull/13) is deliberately separate from this now completed 1,024-photo source.

The earlier tier06 real 8,942-photo source (3,597 true ROI9, 156 strict site×years) is unchanged, the original historical-main ALAN five-year coefficients are unchanged, and MOIS remains the moisture-only repository.
