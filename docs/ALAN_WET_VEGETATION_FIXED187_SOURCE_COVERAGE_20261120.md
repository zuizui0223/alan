# Wetland vegetation dynamics: original matched-pair source coverage

Date 2026-11-20. Supplementary falsification analysis, not independent confirmation.

The existing same-187-pair light slopes are -3.475 without annual open-water trends and -3.498 with them (95% block bootstraps both include zero). Baseline open-water standardized differences +0.722 and +0.745 exceed the preregistered 0.25 threshold; SOURCE_SUPPORT_STOP applies regardless of coefficient stability.

This addition uses repository-pinned original 2000 CORINE 411/412 inland-marsh/peatland points, excluding the coastal-influenced point, and Terra MOD13Q1 250 m 16-day QA=0 NDVI from 2000-2023. An annual NDVI requires data in all three growing-season two-month windows. A site trend requires at least 12 suitable years spanning 18 years, at least two early years (2000-2004), and at least two late years (2019-2023). NDVI represents vegetation greenness, not wetland vegetated area or hydrology.

The 187 original matched pair identities and their 30 blocks remain fixed. The script audits how many of those pairs have vegetation at both arms; NO zero-filling, re-matching, or outcome-based site selection is allowed. Only if at least 30 same pairs in at least 10 blocks have valid NDVI does it estimate nested light coefficients using identical subsamples (baseline, +water change, +NDVI baseline and slope), with spatial-block bootstrap. Failing coverage returns SOURCE_SUPPORT_STOP and no fabricated coefficient.

Crucial limitation: existing wetland NDVI covers only 18 original mapped wetland coordinates, 17 after excluding one coastal-influenced coordinate, versus the 1707 candidate botanical sites from which 187 pairs were built. Full vegetation-change validation across all 187 pairs requires recovering the original all-site location mapping from the unpublished local gbif_annual/annual.duckdb geo table and acquiring CORINE 2000 versus 2018 wetland cover changes or a comparable annual land-cover series for all matched sites. The published 187-pair trends table lacks that complete location mapping; therefore a 187-pair vegetation-adjusted estimate must not be claimed until verified.

Run script: scripts/analysis/alan_wet_vegetation_identical_pairs_20261120.py
Run tests: tests/test_alan_wet_vegetation_identical_pairs_20261120.py


## Executed original-source result (2026-11-20)

The actual source-pinned GitHub Actions run https://github.com/zuizui0223/alan/actions/runs/44666908447 PASSED all 3 new guards and data processing; its immutable result artifact is https://github.com/zuizui0223/alan/actions/runs/44666908447/artifacts/14597874040 (ZIP SHA256 afa3fb2c03e225c7d3c00be3ec80d6f3171b92156af059de58620a47ab7bd489).

- Original 187 development-matched pairs / 30 blocks unchanged.
- 17 long-running QA-valid inland-wetland vegetation sites.
- Within the original 187 pairs, just **1 bright-arm site and 0 comparison-arm sites** have these independent long-running mapped-wetland NDVI observations.
- **0 original matched pairs, 0 blocks** have both wetland vegetation histories. Result: VEGETATION_SOURCE_SUPPORT_STOP; **no vegetation-adjusted ALAN effect estimable**. Do not describe this as evidence ruling out vegetated-wetland loss.
- The previously verified same-187-pair open-water sensitivity still has coefficient -3.475 versus -3.498, confidence intervals spanning zero and failed initial-water balance; no added primary causal claim.

Numerical source-signed output committed at results/all_species_forest/wet_falsification/wetland_vegetation_fixed187_20261120.json.
