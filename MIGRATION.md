# MOIS → ALAN transfer record (2026-11-08)

- **Source:** `zuizui0223/mois`, draft PR [#199](https://github.com/zuizui0223/mois/pull/199), immutable source head commit `63458d81836822126851489352ad29597a87e4d2`.
- **Destination:** `zuizui0223/alan`.
- **Transfer policy:** Copy original research text, scripts, tests, aggregate results, and acquisition manifests without changing the scientific results. Preserve original relative Python module locations. Store old Actions YAML under `archive/mois_workflows/` so bulk data download/extraction cannot accidentally re-run when the new repository receives commits.
- **Original source data:** MOIS GitHub Actions artifacts/releases are the traceable source; the transferred aggregate CSVs/JSONs are *not* the full original images, canonical source observations, VIIRS rasters, or WorldPop products.
- **Scientific baseline:** same-year 2019/2021/2023/2024/2025 island×ALAN photographed white ROI spatial interaction −2.497 pp [−3.949,−1.046] in 7,982 photo cell-years; expanded-cohort −1.205 pp [−2.487,+0.077]. Preserve both estimates and their differing target samples.
- **Reproducibility gate:** Legacy workflows and scripts require source artifacts in `zuizui0223/mois`. Re-run them in ALAN only after validating their inputs, SHA constraints, imports and GitHub Actions permissions; no assertion of a new successful run in ALAN is made by this migration.
- **MOIS scope:** continue global moisture-related photographic colour work there. This transfer does not delete scientific history or invalidate earlier MOIS ecological results.
