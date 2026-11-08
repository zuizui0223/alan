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
