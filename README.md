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
