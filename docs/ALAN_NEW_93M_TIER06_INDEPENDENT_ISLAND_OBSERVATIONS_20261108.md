# ALAN — genuinely new island original observation IDs from the 93,209,351-photo upstream source (2026-11-08)

**SHA-verified source-only acquisition completed:** [ALAN original tier06 source run #42446707775](https://github.com/zuizui0223/alan/actions/runs/42446707775) **SUCCESS**. The original-photo inference step has a separate [run #42447120685](https://github.com/zuizui0223/alan/actions/runs/42447120685) and is **not** part of this data-only source certificate.

## New data, not duplicate old photo angles

The prior 1,181,146 original-photo union contained 148,264 additional photo IDs, but *all* those photos belonged to existing ALAN canonical observation UUIDs, so they could rescue previously unclassifiable observations but add no independent field observations. This source step selects genuinely novel UUIDs:

| Proven authentic source opportunity | Observed count |
| --- | ---: |
| Entire globally indexed original-photo frame | **93,209,351 photo IDs** |
| **Only the immutable tier06**: 16 SHA-pinned eight-part photo-metadata groups (not pixel images) | **597,070 photo IDs** |
| Distinct original observation UUIDs *absent from existing ALAN 1,032,882-observation canonical source*, after selecting one photo per observation | **248,071 UUIDs** |
| New original observations with supported photo year and usable point | **148,516** |
| New observations with actual GSHHG island or mainland polygon assigned | **144,261** |
| New **island** original observations | **17,223** |
| Distinct GSHHG islands represented by novel observations | **861** |
| New island observations with actual WCVP accepted taxon | **17,073** |
| New island observations with actual original-year VIIRS and WorldPop matching the frozen source site | **8,942** |
| Distinct island photo cell-years with these new observations and real original-year matched environment | **4,776** |
| Distinct islands with independent new environment-matched photo opportunities | **369** |

This uses a single pre-published observation-hash tier of the *entire 93m candidate frame*, avoiding rescanning 93m photos or downloading 93m RGB images. It excludes the entire previous 1,032,882 canonical **observation UUID** set *before any new photo selection*. Photo choice within a genuinely new UUID is deterministic and **blind to hue, technical ROI success or flowering-state outcome**. No new colour label is made by this source-only stage.

## First fixed original-photo study

**256 new original photo IDs, 256 distinct observation UUIDs, 256 GSHHG islands and 256 original 0.05-degree photo-cell × year sites** have been SHA-frozen as the first true-photo acquisition queue. Each was observed in a target year (2019/2021/2023/2024/2025), has a genuinely WCVP accepted nominal taxon, a supported original photo URL/format, and matches a real earlier source grid with same-year VIIRS and population. The source queue SHA256 is `7fb7bbb60c556b3a2a31e10bda60715996e6738a662896c377f5f1569ff3bc56`; it is stored inside the source Actions artifact `alan-93m-tier06-256-new-island-observation-original-photo-queue`.

Priority is **island coverage**: at most eight new original observations per island and earlier source-undersampled islands ahead of previously well sampled ones. The first 256 happened to represent 256 different islands; this is an **island-targeted conditional subsample**, not an unbiased worldwide flower-colour incidence survey. The original tier06 observation-hash selection probability `0.0064` does **not** equal the final inclusion probability after the island-first prioritization.

### Scientific stage boundary

The first 256 are, at this stage, **not yet proven flower photos** and do not contribute 256 known flower-white labels to the old strict9 710 island photo-cell-year denominator. Only the [independently executing original-RGB YOLO11/EfficientSAM workflow](https://github.com/zuizui0223/alan/actions/runs/42447120685) can establish how many photo images were actually obtained, which real FCP flower-ROI9 fractions were measured, and which were technical unknowns. The original frozen historical 7,982/710 strict 5-year model remains untouched. A new-source island-first acquisition arm must be distinguished from historical-main in any subsequent ecological fit.

### Provenance

- [Signed source and photo queue](https://github.com/zuizui0223/alan/actions/runs/42446707775): original full93m tier source SHA256 `4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b`, all 16 original tier06 group artifact SHA values pinned and verified.
- [New accepted-taxon / island / year source summary JSON](../results/alan_new93m_tier06_independent_island_photo_source_20261108.json).
- [Exact scripted selection](../scripts/acquisition/alan_new_93m_tier06_independent_island_source_20261108.py), [synthetic source-UUID exclusion tests](../tests/test_alan_new_93m_tier06_independent_island_source_20261108.py).
- Source tiers originate from [the published MOIS full 93m photo index](https://github.com/zuizui0223/mois/releases/tag/mois-complete-93209351-original-photo-target-tiers-v1); **all future acquisition and ecological interpretation remain in ALAN**. MOIS is reserved for moisture × floral colour.

No biological loss, pollinator shortage, native-flora presence or ALAN causal effect is established by a metadata-only source extension.

## All 8,942 now scheduled into 35 SHA-frozen disjoint 256-photo waves

[ALAN source run #42447983019](https://github.com/zuizui0223/alan/actions/runs/42447983019) **SUCCESS**. The already signed 8,942 observations from all 369 source-supported islands were scheduled without any additional photo inference. The exact original first256 sequence and old source image IDs are unchanged, so running FCP inferences and future queues will not accidentally schedule these observations twice.

- **34×256 + 238 = 8,942** unique original photo IDs and unique original observation UUIDs in **35** immutable prospective batches.
- Batch 1: **256 observations in 256 islands** (identical to the first published and currently processing pixel source).
- Batch 2: **256 new observation UUIDs over 249 islands**, favouring islands not previously selected in the first batch. As the available pool is depleted, later batches include repeated *islands* but never repeated observation UUIDs.
- Input original 8,942 source SHA256 `d0b7f7d9f117036d05fbc7ff8b1e62847ef05e7bbc07efc659c9f6cf895d30f9`; input first256 original SHA256 `7fb7bbb60c556b3a2a31e10bda60715996e6738a662896c377f5f1569ff3bc56`.
- Final 8,942-photo schedule SHA256 `e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd`; per-batch source support CSV SHA256 `b1c60c593690faf0f66a350e0eaa64bc375357f16724abf4ef7d8e91d83cbc42`. The workflows use gzip `mtime=0` for byte-level reproducibility.
- The full schedule and island-batch support are available as the `alan-source-8942-novel-island-UUID-disjoint256-batches` artifact (run #42447983019, artifact ID **13604800899**, temporary Actions retention).
- The targeted island-balanced scheduling procedure was frozen prior to pixel results; it is not a global random sample or a substitute for known design weights.

The 8,942 rows are **pre-colour source-only observation targets**, *not* 8,942 true flowering plants or 8,942 measured strict9 ROIs. Keep the earlier 7,982/710 source-photo regression frozen until separately estimated on a scientifically supported expanded cohort.

## Additional field-support observation from the signed source (before pixel inference)

The SHA-verified 8,942-row original-photo acquisition schedule and the signed first256 were also checked against the previously executed ALAN *strict five-photo/three-accepted-species* **island-level photo cell-year support count**. This is only a sample-coverage indicator (not island occupancy by plants):

| Frozen source layer | Islands with **zero previous strict-qualified photo cell-years** | Islands with ≥1 previous strict-qualified cell-year | Genuine new source photo observation UUIDs |
| --- | ---: | ---: | ---: |
| All 8,942 environment-matched original photo opportunities | **281** | 88 | 8,942 |
| First 256 source targets | **256** | 0 | 256 |

Thus the **first 256 colour-blind original photographs target 256 distinct islands with no previously strict-qualified ALAN photographic cells**. They are especially valuable for expanding *geographic source opportunity*, but one photograph from one island does not by itself satisfy the original strict 5-photo/3-species cell acceptance rules. Neither the presence of a photograph nor a WCVP name in the source establishes a verified blooming native flower. Only successful authentic pixel ROI plus sufficient actual species/photographic opportunities would eventually add a statistically comparable cell.

The first source target years are 2019: 27, 2021: 40, 2023: 58, 2024: 57, 2025: 74; first256 hemispheres are 228 Northern and 28 Southern. The master 8,942-photo opportunity pool contains 7,946 Northern and 996 Southern records, so Southern source-photo opportunities remain less abundant and no unweighted hemisphere or worldwide colour prevalence claim follows.

The signed 35-batch schedule's second 256-photo source group covers 249 islands; it never repeats a photo or observation UUID from the first group. All source support counts refer to already observed original ALAN photographic source eligibility, **not** plant extinction or ecological species richness.

## First source-only 256 now actually measured: 92 authenticated ROI9 successes

The independently executed [ALAN original-RGB/FCP frozen ROI-v4 Actions run #42447120685](https://github.com/zuizui0223/alan/actions/runs/42447120685) finished **SUCCESS**. Original 2026-09-27 photo URLs were fetched from iNaturalist's published open-data S3 exactly once for this new work unit, with the unchanged FCP YOLO11/EfficientSAM source and frozen input photo identity, and image bytes were not retained.

| Genuine newly processed original photographs | Count |
| --- | ---: |
| Source images retrieved, original observations, and original distinct GSHHG islands | **256 / 256 / 256** |
| Four-hue classified, with complete original nine-component floral-region ROI fractions | **92** |
| Original ROI9 nonclassifiable due to ROI/flip gate | **156** |
| Original ROI9 uncertain palette composition | **8** |
| Original image retrieval/acquisition failure | **0** |

The original nine-component ROI success is **92/256 = 35.94%**. Unknowns (164) remain unknown, not "nonwhite" or assumed absent flowers. Each of the 92 successful photos is a genuinely new, previously unrepresented canonical **observation UUID**; however every new island was photographed only once in this first batch, so **none of the 92 is automatically an eligible new strict 5-photo/3-accepted-species island photo-cell**. The frozen historical-main original five-year model and ecological effect sizes remain unchanged.

The signed original colour data artifact in GitHub Actions is `alan-256-genuine-new-independent-island-original-image-ROI9` (artifact ID **13604781981**); its source-photo nine-palette CSV has SHA256 `7789c5aa1916e5b1ff915d217744a80b990795cb24e362cff20007dda96ee94b`. Method fingerprint: frozen FCP source Git blob `b1d7091db44682559bbe4b48a9f3b8db648501a6`, unchanged original model source commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`. Numerical receipts and aggregate provenance are also in [ALAN first256 signed results](../results/alan_first256_new_island_real_original_roi9_20261108.json).

**Next acquisition axis:** progress through the already source-authenticated 35 disjoint 256-photo batches, retaining the first 256 image outcomes instead of downloading those images again. Each additional batch must pass true photo/observation-ID uniqueness, image-byte ephemerality, frozen model fingerprint and unknown-colour quality gates. Because the queue is prioritized by island/geographic source coverage, observational source selection differs from historical-main, and no unweighted merged worldwide biological colour rate is estimated.
