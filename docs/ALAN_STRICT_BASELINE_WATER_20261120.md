# ALAN matched wet plants: original-water baseline common-support re-test

2026-11-20 | Pre-outcome falsification design, not independent confirmation.

Scientific problem: the original 187 pairs/30 blocks have open-water baseline
standardized differences +0.722 area and +0.745 water detection frequency.
These violate the fixed 0.25 standardized-balance threshold. Hence the same-pair
slope -3.475 -> -3.498 with annual water-trend covariates is only posthoc
sensitivity, not an identified ALAN-specific causal or mediation effect.

Additional test using the EXACT source SHA-verified 1707-site original frozen
1999-2023 data, never ecological outcomes in matching:
- Repeat the same Hungarian within-100-km-block, no-replacement match.
- Require light-change difference >=1 DN/decade, 0.5 IQR built-change
  caliper, 1 IQR baseline built/light calipers, survey start <=5y,
  timespan <=10y, and the same unchanged standardized distance.
- The PRIMARY new pre-exposure restriction is baseline water-area and
  detection-frequency absolute differences each <=0.5 x scale, and stricter
  <=0.25 x scale is prespecified SENSITIVITY (not tuned on effects).
  The scale is max(empirical IQR, 0.1 empirical SD); this handles
  strongly zero-inflated baseline water, for which raw IQR can be zero.
  Each fraction multiplies the SAME robust scale so the two calipers
  cannot collapse to identical thresholds due to the zero-IQR floor.
- Independently audit the effective number of pairs and spatial blocks,
  actual SMD on built, light0, initial water, climate, and original guard
  (>30 pairs, >10 spatial blocks, every critical |SMD|<=0.25).
- For each fixed within-variant matched cohort, compare the light-difference
  coefficient with/without ANNUAL mapped surface-water trend in identical
  pairs, with 100-km-block bootstrap intervals.
- Evaluate new pair-selection effects as a different subpopulation; the
  new slopes MUST NOT be compared to 187-pair coefficients as an estimate of
  how much water mediates night illumination.

Source is a postselected 58 wet-decline candidates versus 305 other wet plants,
so any effect is exploratory and may be affected by outcome selection and
nonrandom long-term botanical survey design.

This test requires no mapping coordinates or new vegetation data. It does not
claim that wetlands were unaffected: vegetated wetland loss, drainage, hydrology
and nighttime spectral exposure remain unobserved. See other main diagnostics
for the 0/187 verified geographic overlap with original 2000-2023 wetland NDVI.


## Executed source-authentic numerical results (2026-11-20)

Four unit tests and the frozen-original-source numerical replay passed:
https://github.com/zuizui0223/alan/actions/runs/44672675032
Signed ZIP SHA256 f0a68917d2249eb59a18488c08c6545ce86780c59f0629d50eecc25f9a2f2fa6.
Full transparent coefficients, matching IDs digests, critical balance, block-bootstrap CIs:
[strict baseline-water report](../results/all_species_forest/wet_falsification/strict_baseline_water_20261120.json).

| Source audit | Original development-matched | Baseline water caliper 0.5 × robust scale | Baseline water caliper 0.25 × robust scale |
|---|---:|---:|---:|
| Matched source pairs / spatial blocks | 187 / 30 | **133 / 26** | **128 / 25** |
| Standardized baseline water area difference | +0.722 | **+0.015** | **+0.010** |
| Standardized baseline water frequency difference | +0.745 | **+0.017** | **+0.008** |
| Paired recording-rate mean (pp / decade) | -2.476 | -1.896 [-4.174, +0.235] | -1.687 [-4.149, +0.337] |
| Conditional slope without annual water trends (pp / 10 DN) | -3.475 | -1.149 | -0.977 |
| Conditional slope with annual water trends (pp / 10 DN) | -3.498 | **-1.186 [-4.323, +1.893]** | **-1.067 [-4.037, +1.968]** |
| Baseline-water source common-support gate | STOP | **EXPLORATORY_COMMON_SUPPORT** | **EXPLORATORY_COMMON_SUPPORT** |

Here the two actual baseline water calipers are distinct: area thresholds 0.827 / 0.413 percentage points; water-detection frequency thresholds 0.848 / 0.424 percentage points. The original source's 0.5/0.25 calipers were initially inadvertently identical due a zero-IQR fallback under the matched eligible sample, which was detected and corrected in tests before freezing these numerical results. Data-table marginal IQRs are nonzero, but the matching eligible subset can be zero-inflated. The two final strict samples are different.

**Interpretation**: Post-outcome-selection data and CIs spanning zero preclude a light-specific confirmed effect. The sign remains negative with more comparable *initial* water. Its point magnitude decreases when different sites are matched, but that between-sample contrast does **not** estimate hydrological mediation or causal attenuation. The balance gate checks baseline covariates; contemporaneous open-water trends remain unbalanced and adding their slopes within a fixed 133/128 cohort barely changes the point coefficient. Crucially, neither exercise controls **vegetated-wetland area** on the original 187 pairs (0 original paired sites with independently measured NDVI), so source verification is still required before claiming that explanation has been excluded.

Methods and quantities are supplemental Results only, not grounds to revise the primary discovery-stage wet-affinity +0.534 candidate trait contrast.
