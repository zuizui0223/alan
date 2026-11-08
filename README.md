# ALAN — artificial light at night × island floral colour

**Research home for the ALAN × island/mainland floral-colour project, migrated from [MOIS](https://github.com/zuizui0223/mois) on 2026-11-08.** MOIS is reserved for global floral colour × moisture.

## Central question

Does the spatial association of artificial light at night (ALAN; VIIRS upward radiance) with the *photographed flower-region white-colour fraction* differ between island and mainland sites?

## Source-verified result carried forward

- Five observation years: **2019, 2021, 2023, 2024, 2025**.
- **7,982** source-backed 0.05° photo-cell × year records, including **710** island cell-years.
- Island × log(1+VIIRS) spatial interaction: **−2.497 pp** (nominal spatial-cluster 95% CI **−3.949 to −1.046**; p **0.000802**), after climate, elevation, population, photography effort, geographic and observation-year controls.
- On the broader observed-photo eligibility cohort (20,543 cell-years), the estimate attenuates to **−1.205 pp** (95% CI **−2.487 to +0.077**); matched same-year/geographic comparisons and repeated sites do not establish a universal island-specific light response.
- The measured endpoint is an **original photographic floral-ROI white component**, NOT independently verified white-flowered species prevalence, causal trait evolution, effective nocturnal pollination, or local extinction.

The original source runs, numerical results and caveats are preserved in [the migrated draft](manuscript/MOIS_ISLAND_ALAN_OBSERVATIONAL_COLOUR_ECOLOGY_20261108_RESULTS_DRAFT.md) and [analysis notes](docs/MOIS_ALAN_ISLAND_FIVEYEAR_EMPIRICAL_20261107.md).

## Structure and execution

- `manuscript/`: draft observational Results/Discussion.
- `docs/`: primary empirical and data-source audit notes.
- `scripts/`: original research scripts, **preserved under their existing import paths**.
- `tests/`: historical tests; successful original GitHub run links remain in the source notes.
- `results/`, `config/`: published aggregate tables, checksums and source queues; **not the original photo/VIIRS data archives**.
- `archive/mois_workflows/`: copied legacy GitHub Actions as **non-executing historical recipes**. They refer to MOIS action artifacts and old branch triggers; do not simply move them into `.github/workflows/` without adapting and validating external inputs.

## Migration / provenance

See [MIGRATION.md](MIGRATION.md). The original signed image/data artifacts remain in the source MOIS release and Actions history; source files have not been remeasured and missing colour values have not been filled. Scientific analysis continuation belongs here, and moisture–flower-colour analysis remains in MOIS.

## Native ALAN reproduction (verified)

On 2026-11-08, ALAN's own [five-year source-verified GitHub Actions run](https://github.com/zuizui0223/alan/actions/runs/42438862639) succeeded using the historical MOIS original-photo and same-year VIIRS/WorldPop source artifacts. It reproduced **7,982 photo cell-years (710 island)** and the island × log(1+VIIRS) white-ROI spatial association **−2.497045 pp**, 95% CI **[−3.948547,−1.045542]**. The stricter same-year same-geographic island/mainland comparison was **not individually resolved**; no light-caused biological white-flower loss is inferred.

Use `.github/workflows/alan-native-fiveyear-replay.yml` via the **Actions → Run workflow** manual control on the default branch. This workflow intentionally leaves the historical MOIS workflows inactive. See [full original-run outputs, safeguards, and archival deadlines](docs/ALAN_NATIVE_REPLAY_VERIFIED_20261108.md).

## Exploratory island-area result (2026-11-08)

A source-verified new ALAN run on real GSHHG area and original flower-region ROIs found a **positive area × nightlight interaction** for within-island photographic whiteness: **+2.209 pp per tenfold island area** (nominal island-cluster 95% CI **+0.007 to +4.411**, 44 spatially supported islands/650 photo cell-years) in the strict cohort, but **+0.998** (CI **−0.390 to +2.387**, 95 islands/1,884 cells) in the broader measured-photo cohort. The source-support audit found that the larger-area island third supplied ~84% of photographed cell-years in both cohorts. Giving **each island equal total weight** produced **+5.152 pp** (island-cluster nominal CI **+2.559 to +7.744**) in the strict 44-island cohort and **+3.385 pp** (**+1.319 to +5.450**) in the broad 95-island cohort. These estimates support a **positive observational area moderation** under both island-equal cohorts, but are exploratory and not independent causal proof. Thus a smaller-island stronger-negative-light slope is **exploratory and sensitive to cohort support**, not a universal ALAN effect or confirmed pollinator mechanism. See [original empirical results](docs/ALAN_GSHHG_ISLAND_AREA_LIGHT_EMPIRICAL_20261108.md) and [successful ALAN run](https://github.com/zuizui0223/alan/actions/runs/42439477013).

## Island area versus mainland distance (verified)

The new source-verified [coastline-isolation analysis](docs/ALAN_MAINLAND_ISOLATION_VS_ISLAND_AREA_EMPIRICAL_20261108.md) compares island size with nearest GSHHG mainland coastline distance, using **the same 2019–2025 original photographed floral ROI source**. In joint within-island, **island-equal** models, the island-area × nightlight interaction remains positive in both the strict **44-island** cohort (**+5.384 pp**, 95% CI **+2.678 to +8.089**) and broad **95-island** cohort (**+2.939 pp**, 95% CI **+0.833 to +5.045**). The mainland-distance interaction is uncertain in both. **Cell-year-equal area intervals also include zero**; this remains a sampling-sensitive, exploratory spatial association, not causal evidence for flower evolution or nocturnal pollination changes. [Reproducible ALAN run](https://github.com/zuizui0223/alan/actions/runs/42440394223) · [signed result summary](results/alan_isolation_area_light_20261108_source_verified_summary.json).

## Photographed taxon composition and island area × nightlight (new real result)

ALAN has now [source-verified and executed a same-island, same-photo-cell-year decomposition](docs/ALAN_ISLAND_AREA_TAXON_COMPONENTS_EMPIRICAL_20261108.md) for the positive island-area × nightlight observation after adjusting continental isolation. In the **equal-island** model, strict 44-island interaction **+5.384 pp = +4.296 pp** from observed accepted taxa's *same-photo-archive* white-reference component **+1.088 pp** from taxon-relative photographed colour; broad 95-island **+2.939 = +2.059 +0.880 pp**. These are exact algebraic photography-based splits, **not biological species replacement or causal mediation**.

A more independent **2019/2021/2023 photographic species-colour reference evaluated in 2024/2025** has support in **29 strict / 50 broad islands**, but the island-equal earlier-year taxon-reference component has intervals overlapping zero. Thus **species assembly cannot yet be identified as the unique ecological mechanism**. Within-taxon photographic deviations remain possible but include observer/season/ROI effects, not verified flower-colour evolution.

[ALAN execution #42442841980](https://github.com/zuizui0223/alan/actions/runs/42442841980) · [Compact source-signed result](results/alan_island_area_taxon_components_20261108_verified_summary.json). The earlier island-area and coast-isolation results remain unchanged.

## Real additional floral ROI information from already measured sibling photos (2026-11-08)

The signed MOIS 1,181,146-original-photo / 451,643-strict9 archive contains **148,264 extra photo IDs** beyond ALAN's canonical 1,032,882 original one-photo-per-observation reference. They are all **other photos of existing observation UUIDs, not new independent observations**. Of 49,940 genuine extra strict9 photographs, **18,950 existing observation UUIDs** have a strict9-successful sibling when their original representative photo was technically unknown.

An ALAN-native SHA/GSHHG/WCVP/source-year audit found that admitting these genuine sibling-photo colour results **without increasing photo-opportunity or observation counts** increases **strict-eligible photo cell-years 7,982 → 8,515 globally and 710 → 793 on islands (+83, spanning 25 islands)**; the 3-photo/2-accepted-species cohort rises **20,543 → 21,680 globally and 2,026 → 2,155 on islands (+129 across 41 islands)**. Forty geographical cells with inconsistent original accepted-taxon reconciliation were excluded from additional admission, not removed from the prior results. **No new island×nightlight regression was fitted on the expanded source:** the values are genuine measured-photo opportunity gains, not a new ecological effect.

[Complete real source evidence](docs/ALAN_MORE_AUTHENTIC_SIBLING_PHOTO_ROI_SOURCE_20261108.md) · [Successful ALAN data-only run](https://github.com/zuizui0223/alan/actions/runs/42445111291) · [aggregate JSON](results/alan_20261108_extra_measured_sibling_ROI9_photo_opportunities.json).

**To add truly independent photo observations next**, select new observation UUIDs from the already-indexed 93.2m-photo upstream candidate frame before downloading any image pixels; keep source-reference and targeted refill selection separate. The MOIS repository remains moisture-focused.

## Real acquisition from the full 93.2m photo index (2026-11-08)

The original-photo source is now expanded beyond ALAN's old one-photo-per-observation canonical cohort. Source-only [run #42446707775](https://github.com/zuizui0223/alan/actions/runs/42446707775) authenticated **597,070** records from one pre-signed original 93.2m photo-index tier, excluded all prior canonical observation UUIDs, and identified **248,071 new observation UUIDs**, including **17,223 on 861 GSHHG islands**. **8,942** new island source observations on **369 islands** already have accepted taxonomy and actual same-year VIIRS/WorldPop source support. Their original photo identities have been frozen into **35 disjoint batches** (34×256 + final 238; [signed batch run](https://github.com/zuizui0223/alan/actions/runs/42447983019)).

The first 256 original photographs on 256 different source-supported islands **have been genuinely acquired and processed** with the unchanged FCP original YOLO11/EfficientSAM ROI-v4: **92 photographs yielded numerical 9-component ROI**; 164 were technically unclassifiable (156 ROI/flip gate, 8 ambiguous palette), and **zero photo downloads failed**. [Original signed image result](results/alan_first256_new_island_real_original_roi9_20261108.json) · [successful source-photo inference](https://github.com/zuizui0223/alan/actions/runs/42447120685). These new photos are distinct observation UUIDs, but one photo/island does not meet the original strict 5-photo/3-species ALAN cell definition. The previous ecological regression remains frozen; the new acquisition is a separately labelled, island-prioritized source cohort, not a globally representative white-flower sample.

## Independent real island photographs — first five acquisition batches

The island-prioritized 93m-source acquisition arm has independently completed original RGB+ROI-v4 inference on **1,280 genuinely new observations**: **492 true numerical original nine-component floral-region ROIs** and **788 technical unknowns**, zero failed original-image downloads, [combined source signature verification](https://github.com/zuizui0223/alan/actions/runs/42456790465). The 1,280 photo/observation UUIDs cover **369 GSHHG islands and 1,060 original photo cell-years**, but just **one** cell has >=5 original photo opportunities and >=3 truly measured WCVP accepted species. This is a different, deliberately island-rich sample from the original global 7,982/710 strict cell-year panel; it is not an increase of 492 regression cells. [Scientific data intake note](docs/ALAN_REAL_FIRST1280_ISLAND_OBSERVATIONS_ROI9_SOURCE_20261108.md) · [aggregate numerical result](results/alan_real1280_original_new_island_ROI9_verified_20261108.json).
