# ALAN 2026-11-15 | Li original plants: the positive light response is baseline-zero sensitive

## Source and execution

The original author [Dongbo Li et al. 2026 bioRxiv preprint, DOI 10.64898/2026.07.10.737511](https://doi.org/10.64898/2026.07.10.737511) [seedset.csv](https://github.com/Dongboli/experimental-data/blob/main/Artificial%20lights%20at%20night/seedset.csv) Git blob `b616aca5beeba1457549841477d61ca4f796a147` was re-downloaded and checked with the Git blob SHA1 algorithm. Original 79 plant rows, 78 genuine known seed outcomes. **The one missing seed measurement is in the point-light P arm**, not silently assigned zero in C or L.

New independent source matching holds field (3 true light experimental meadows) × plant species (ML/SN/ZC) × field flowering round constant, keeping only stratum pairs where both L (linear array of 26 lamps) and C (unlit control) had true original measurements. This yields **12 matched source strata, NOT 12 independent light allocations**. The source-author light geometry/number confounding is unchanged.

The [exact-source Actions run #43720758472](https://github.com/zuizui0223/alan/actions/runs/43720758472) **SUCCESS, 5/5 local tests**. Results persist in [source-signed compact summary](../results/alan_Li2026_matched_baseline_floor_source_sensitivity_20261115.json) and [original-source script](../scripts/analysis/alan_Li2026_matched_field_taxon_round_20261115.py).

## Exact new findings

| Source population subset (post hoc) | Matched strata | Difference in seeds/initial bud, L−C, equal strata | L/C seeds per starting bud, matched pooled | Strata sign (+/−) | True field mean sign (+/−) |
|---|---:|---:|---:|---:|---:|
| All identical field × species × round | **12** | **+5.805** | **3.608×** | **10/2** | **3/0** |
| Require ≥2 original plants in BOTH treatment arms | **9** | **+6.675** | **3.633×** | **7/2** | **3/0** |
| Dark C had measured seed count >0 | **6** | **+6.364** | **2.790×** | **4/2** | **2/1** |
| Both ≥2 plants /arm and dark C seeds >0 | **5** | **+6.917** | **2.819×** | **3/2** | **2/1** |

**Key newly demonstrated sensitivity:** **6/12 source-matched strata have zero observed seeds in the dark control, and all 6 had a positive L−C contrast.** Once the control-zero strata are excluded, field 2 changes to a *negative* mean seed/initial-bud difference (**-1.503**) and field signs are 2 positive / 1 negative, rather than 3/0. Requiring ≥2 plants per arm leaves nine strata and three positive field means, but pairing BOTH stricter conditions leaves only five strata, of which three are positive, two negative. Hence the original whole-source and all-strata signals are truly positive descriptive associations, **but 'positive under lighting across all environments' is not supported by this baseline-sensitive analysis.**

The pooled ratio decrease from previously reported all-source 4.239 to matched 3.608, and to 2.790 with dark-zero exclusion, is **not** a controlled treatment adjustment. The dark-zero exclusion explicitly selects groups by their *observed outcome* and should **never** be presented as causal unbiased sensitivity, a policy recommendation to illuminate dark habitats, or evidence that a light source rescues pollen-starved plants.

**Important distinction:** fruit per starting bud was still positive in all three field-equal averages even after dark-zero exclusion (and in 5/6 of those source strata). But fruit production, seed-positive individual, seed per fruit, and seed per original open flower are different estimands. The possibility of a 'reproductive establishment hurdle' is a mechanistic HYPOTHESIS, not proof of moth pollination, corridor effect or autonomous-selfing substitution.

## Consequences for 9-paper ALAN evidence synthesis

ALAN's reproductive association can have stage-specific sign reversals and **zero-inflated / hurdle-like response** among species, fields, light configurations and flowering rounds. The correct meta-analysis design should predefine **at least three parallel endpoints**: (1) probability that an originally assigned flower/plant produces any fruit or seed, (2) total seed number per original opportunity (seeds/original bud/flower/plant, clear denominator), and (3) seed quality conditional on formed fruit (distinct, post-treatment selected). Include actual randomized light site, year, and plant ID to avoid treating 12 or 79 rows as lighting replicates. Sample species and lighting spectra/configurations without selecting only negative studies.

There remain **nine** source-confirmed directly manipulated plant-reproductive ALAN publications, plus three observational gradients; this is ONE of nine, not 12 new studies. The global `β_V` of ALAN reproduction by independent same-population marginal nighttime reliance remains **UNESTIMATED (k=0 fully prepared)**. New variation in source endpoints **does not** imply a usable meta regression or island × ALAN experimental reproduction result. Future research should prioritize original lit/dark × day/open factorial cells (Lu 2026, Macgregor 2019), then rederive cluster-aware covariances.

## Provenance and inferential controls

Data source is original published-author CSV with Git blob identity, not visual-only or simulated values. Source tests verify source ID and original treatment identities; exclusions are displayed, no imputed seeds. A 3-field unanimous sign has a mathematical two-sided sign probability **0.25** under simplistic independent exchangeable field signs, **not a verified randomization p-value**. The source does not certify equal light wattage across linear/point configurations. No raw individual author rows are redistributed.
