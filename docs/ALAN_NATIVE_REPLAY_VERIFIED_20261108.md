# ALAN native execution — independently re-run in the new repository (2026-11-08)

## Execution outcome

- **New ALAN Actions run:** [42438862639](https://github.com/zuizui0223/alan/actions/runs/42438862639), **SUCCESS** (all steps); [uploaded numeric source-receipt artifact](https://github.com/zuizui0223/alan/actions/runs/42438862639).
- Native workflow: `.github/workflows/alan-native-fiveyear-replay.yml`.
- Source code: copied byte-for-byte from frozen MOIS [PR #199](https://github.com/zuizui0223/mois/pull/199) head `63458d81836822126851489352ad29597a87e4d2`.
- The original MOIS GitHub Actions source artifacts were **actually downloaded from ALAN** using the ALAN workflow token. Source ZIPs and 2019/2021/2023/2024/2025 output receipts were checked; the 2025 independent known-land source SHA matched the original source, with the previously documented benign full auxiliary geometry-area rounding difference.
- Original images downloaded or reinferred: **0**. Technical unknown ROI colours imputed: **0**.

## Numerical reproduction in ALAN

| Result | ALAN value |
| --- | --- |
| Real photo cell × year observations | **7,982** |
| Island cell × year observations | **710** |
| Geographic blocks (clustered over years) | **317** |
| Islands with >=1 original qualifying cell-year | **101** |
| Islands with enough spatial variation for within-island FE | **44** |
| Global island × same-year log(1+VIIRS) slope difference | **−2.497044599606764 percentage points** |
| Nominal geographically clustered 95% CI | **[−3.9485472403480117, −1.0455419588655164]** |
| Nominal exploratory p | **0.0008023343965928109** |
| Same-year, same 5° block island/mainland subset (1,553 observations, 119 islands) | **−2.6639915045756295 [−6.283405991042761, +0.9554229818915014]**, p=.14149 |
| Within-island FE (44 islands) | **−1.7054134634487526 [−3.953793521494526, +0.5429665945970206]**, p=.13342 |
| Largest island excluded | **−2.709248479669408 [−4.317315927800733, −1.1011810315380832]**, p=.001024 |

Every one-year-out point estimate was negative (2019 −2.574; 2021 −2.141; 2023 −2.370; 2024 −2.246; 2025 −3.243). These are overlapping photographic subsamples, not five independent biological replications.

**Interpretation:** ALAN now executes the existing five-year spatial association independently of MOIS code checkout. This does **not** establish white-flowered plant losses, causal nighttime-light selection, effective moth pollination or within-species colour evolution. The stricter same-year same-geography comparison and within-island effect both have confidence intervals crossing zero. The response is original photographic nine-palette *flower-region* white mass, averaged equally by photographed accepted taxa within source photo cells.

## Re-execution and retention

Once the workflow has been merged into ALAN's default branch, run it via **Actions → ALAN native original-photo five-year island-nightlight replay → Run workflow** (manual `workflow_dispatch`). It is deliberately **not configured to auto-run the expensive source replay on every main-branch update**. Its push trigger serves the temporary test branch only.

The earlier GitHub Actions files remain **inactive under `archive/mois_workflows/`**; never move them back into active workflow paths without checking permissions and source. The ALAN execution still depends on historical MOIS artifacts rather than maintaining a separate copy of full original records. The earliest critical historical source expires **2027-01-17**, later photo and raster artifacts in January–February 2027; the new ALAN outcome artifact expires **2027-02-06**. For permanence, source archives would need a separately approved durable archiving strategy that preserves source rights and checksums.

The build/result metadata are reproducible; p values above are nominal from exploratory models and should not be promoted as new preregistered confirmation.
