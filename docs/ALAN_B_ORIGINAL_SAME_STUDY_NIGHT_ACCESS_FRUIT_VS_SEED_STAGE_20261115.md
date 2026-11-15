# ALAN B: Paired source diel fruit-vs-seed stage dependence (2026-11-15)

## Reproduced original-data result

The historical source file results/night_exclusion/experiments.csv contains 42 rows with daytime-only versus day+night-access reproductive comparisons. Original pairs were matched by study, taxon, locality, site period, year, cultivar, and additional treatment. **14 same-stratum pairs** have both a fruit endpoint and a seed endpoint.

Four pairs were removed from the conservative analysis for documented source-denominator reasons, not for effect sign: the fruit source values in Anderson_2015 (two years) and Scopece_2018 exceeded 1, so are not certified as fruit set probabilities; the 2025 French small-fruit quality table does not match formed fruit counts in 9 of 20 taxon-year-treatment cells, so the Rubus stage pair cannot be treated as safely aligned. Both the broad 14 and strict 10 are preserved.

**10 source-quality stage pairs from 6 original papers remain**; all 10 have larger marginal night-access loss for the recorded seed endpoint than the fruit outcome, V_seed > V_fruit where V=1−D/O. This is descriptive, and none is a new ALAN lit/dark plant reproductive experiment.

| Study | Plant | Strict site/year/sex strata | Mean across strata: V_seed - V_fruit |
|---|---|---:|---:|
| Stone 2018 | Silene vulgaris | 2 | **+0.5162** |
| Morse 1983 | Asclepias syriaca | 4 | **+0.4220** |
| Robertson 2021 | Malus domestica | 1 | **+0.1295** |
| Jabbari 2013 | Saponaria officinalis | 1 | **+0.1254** |
| Vanderplank 2020 | Silene nutans | 1 | **+0.0704** |
| Noone 2022 | Vaccinium angustifolium | 1 | **+0.0362** |

Original fruit/seed effects within studies are correlated; 10 sites are not 10 independent studies, and a seed/formed-fruit outcome is conditional on fruit formation, sometimes using a different denominator from fruit set. The exact original sampling variances, plant-level covariance and any missingness selection are not in this legacy extracted table. There is **no valid confirmatory confidence interval/p-value** for the six source means from this table, nor proof that larger seed-stage effect means a stronger total effect on initial-flower seed production or survival.

The most useful ecological insight is that nighttime access dependence is **reproductive-stage-specific**; a binary fruit formation endpoint may appear weakly affected even when seed traits are more responsive. Conversely, seed quality among surviving fruits is post-treatment selected, so cannot be multiplied into unobserved total seed yield across all originally assigned flowers without source links.

### Reproducibility and boundaries

[GitHub Actions #43706854701 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43706854701), six source integrity tests; artifact #14159145500, ZIP SHA256 9a69ea78ddec4acfe7902e07df59e83e99d71150b0502d8d2b5a977a91873fa2.

Code: scripts/analysis/alan_B_original_DvN_stage_contrasts_20261115.py, tests/test_alan_B_original_DvN_stage_contrasts_20261115.py.
Output: results/night_exclusion/stage_marginal_diel_pairs_20261115.json and stage_marginal_strict_pairs_20261115.csv. Historical source unchanged.

This is B source evidence for the pathway composition and endpoint choice, not the desired A effect: no same-population artificial lighting contrast, no independence-certified V and y covariance. Confirmatory ALAN meta regression beta_V and island-specific ALAN reproductive interaction remain **not estimated**.
