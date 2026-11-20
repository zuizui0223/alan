# ALAN wet-habitat: distinguish illumination from open-water change and settlement growth

2026-11-19 | exploratory falsification extension | **not a new independent positive result**

## Why this analysis exists

The main species-composition result (+0.534 moisture difference in 158 matched candidates versus 1,419 references) can arise from coincident wet-habitat decline, urbanization, survey selection, and observer effects; it is not proof that night illumination causes loss of wet plants. The stricter NPMS species × habitat × year test is sensitive to observation-status restrictions. This file records two falsification paths and a pre-fixed *common-support* test.

## Test 1: direct surface-water change

Existing ALAN source work already measures GLAD Landsat open-water cover and water-detection frequency around **1,707 long-term botanical survey coordinates (1999–2023)**. It compares observed declines between candidate wet plants and other wet plants, conditional on temperature, rainfall, building cover, field effort, and nearby open-water history. In the 425 open-water-"maintained" coordinates, 108 brightening coordinates gave a candidate minus comparison trend difference of **−0.85 percentage points per decade [−2.41,+0.70]**. Spatial-block held-out random-forest improvements from adding light *after* water were −0.95% and +0.36% of baseline MSE in two splits. These are not corroboration of a light-specific driver.

A second analysis of **17 CORINE-mapped inland wetland/peatland coordinates** measured annual GLAD water transitions and NDVI; only **one** was brightening over the full period. It could not identify water-adjusted ALAN effects in wetland sites.

- GLAD Water Dynamics v2: https://storage.googleapis.com/earthenginepartners-hansen/waterC2/download.html (1999–2025 as currently documented).
- JRC GSW YearlyHistory v1.4: https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_YearlyHistory (1984–2021, 30 m, seasonal and permanent water classes).
- Under-canopy inundation, soil moisture, drainage ditches, and wet grassland loss are not equivalent to open-water extent. Treat no mapped surface water as *unmeasured hydrology*, not stable wetland.
- A JRC replication must use **1999–2021 intersection** and the same coordinate windows and water-quality filters, then compare to GLAD; do not splice 2022–2024 with nonexistent JRC data. Pre- and post-change windows must precede the botanical response period to avoid conditioning on downstream changes where possible.

## Test 2: locations with similar settlement growth but different light trends

A previously frozen 1999–2023 paired wet-plant site table was created in commit 6a9a484:
\`results/all_species_forest/wet_water/coordinate_trends.csv.gz\`.
It contains site-specific annual-nightlight trend \`light\`, GHSL-derived built-up change \`built\`, baseline built-up \`built0\`, baseline light \`light0\`, water cover and detection change, survey start/span, climate and observer-effort trends, and **candidate-minus-reference recording-rate trend** \`rate_difference\`.

The script \`scripts/analysis/alan_light_development_overlap_20261119.py\` works without downloadable private source caches: it fetches the *fixed historical Git blob* through \`git show\` and verifies its SHA256. Its matching procedure was written before reading any outcome:

1. Require light slope >0 in the brighter arm and <=0 in the comparison arm, with a difference >=1 DN per decade.
2. Match within the same 100-km block, without replacement; require built-up change within 0.5 observed IQR, baseline built and baseline light within 1 IQR, survey start within five years, and span within ten.
3. Among eligible pairs, minimize standardized distance across development, baseline light, survey time, climate, observer-effort trends and baseline water. Additional water sensitivity also matches annual water area and detection-frequency trends.
4. Report actual number of matched pairs and independent 100-km blocks, standardized built and light imbalances, paired recording-rate trend differences, and 100-km-block bootstrap intervals. If fewer than 30 pairs, 10 blocks, or critical covariate SMDs exceed .25, set **SOURCE_SUPPORT_STOP**. The code still prints descriptive estimates but these must not be treated as evidence when stopped.

**Important limitations**: these botanical outcomes compare preselected 58 candidate wet plants with 305 other wet plants **in the same discovery-era records**; they cannot validate that selection independently, nor reconstruct all 3,901 species. The GHSL settlement change is based on a 1990–2020 epoch contrast rather than a fully contemporaneous annual panel for 1999–2023. Upward satellite radiance measures neither ground illumination at a plant nor light spectrum. Built-up change is not the same as drainage/management. The within-block matched difference can reflect unmeasured differences in survey protocol, land management, selection and plant community composition. This is a **common-support and falsification diagnostic**, not an ALAN treatment effect.

## Prespecified decision

- If matched pairs are few or unmatched on development: **no claim of separable nightlight effect**.
- If matching works but candidate-minus-reference light contrast includes zero or flips after water adjustment: **downgrade the wet-species observation to an indicator of coincident landscape change**.
- If a contrast persists after both forms of adjustment: treat it as **exploratory until an independent botanical dataset or held-out epoch replicates it**.
- Do not merge into \`main\` as the primary causal mechanism; PR review stays draft until source results and independent validation are assessed.

## Executed source-verified audit and outcome-independent STOP

[Successful GitHub Actions run #44528843759](https://github.com/zuizui0223/alan/actions/runs/44528843759) recomputed the frozen historical source (raw gzip SHA256 \`4d151698e0da5ee754326220d3bb817b141da00770299e7f8910366d8771e74f\`). All source tests passed. Full audit: [JSON](../results/all_species_forest/wet_falsification/light_vs_development_overlap_20261119.json).

| Source comparison | Development-matched | Development + water-change-matched |
|---|---:|---:|
| Pairs / independent spatial blocks | **187 / 30** | **114 / 25** |
| Median matched light-slope contrast (DN per decade) | 4.58 | 5.80 |
| Mean site-paired candidate-minus-reference recording-rate trend difference (pp per decade) | −2.476 | −1.020 |
| Source-block bootstrap 95% descriptive interval | [−5.291,−0.237] | [−2.987,+1.653] |
| Baseline open-water area standardized mean difference | **+0.722** | **+0.360** |
| Baseline water-detection frequency standardized mean difference | **+0.745** | **+0.360** |
| Built-change standardized mean difference | +0.033 | +0.023 |
| Source-support gate | **SOURCE_SUPPORT_STOP** | **SOURCE_SUPPORT_STOP** |

The initial development-matched negative mean is *not a qualified ALAN-specific estimate*: baseline water is poorly balanced, and the reference plant group already belongs to a previously selected candidate survey panel. Adding water-change matching changes the membership from 187 to 114 pairs and attenuates the mean. That comparison cannot attribute attenuation to surface water because the underlying matching population also changes.

### Identical-pair sensitivity without reselecting sites

On exactly the same 187 development-matched pairs, retaining baseline water, climate, light and built controls in both regressions:

- Without annual water-change slopes: relative recording-rate trend slope **−3.475 pp per 10 DN per decade**, block-bootstrap 95% descriptive interval **[−7.287, +0.598]**.
- With annual water-area and frequency change slopes: **−3.498**, interval **[−7.052,+0.610]**.
- Difference −0.023; the water-change covariates **did not attenuate** the matched-pair *slope* on this same sample.

This was added as a post-hoc sensitivity to distinguish change in model covariates from change in matched pair composition, not a pre-registered confirmatory estimator. It has wide bootstrap intervals and the development-matched sample still fails the initial-water balance gate. No evidence here establishes that light, rather than drainage or correlated land-cover change, caused wet-plant decline. Conversely, failure of an open-water proxy to change the point estimate cannot exclude subsurface or vegetated-wetland loss.

**Disposition**: empirical matching was executed, but both matched analyses are **held out of the manuscript's principal evidence**. Do not interpret the earlier negative descriptive interval as an ALAN-specific ecological mechanism. An independent 3,901-species continuous-response model or held-out botanical survey paired with direct hydrology and contemporaneous urban/light trends is the next useful test.

## Reproducibility

\`\`\`bash
python -m pip install numpy==2.1.3 pandas==2.2.3 scipy==1.15.2
python -m unittest discover -s tests -p test_alan_light_development_overlap_20261119.py
python scripts/analysis/alan_light_development_overlap_20261119.py --out /tmp/alan_light_development_overlap.json
\`\`\`

Source context: \`results/all_species_forest/wet_water/README.md\` at commit 6a9a484; \`results/all_species_forest/wet_annual/README.md\` and \`results/all_species_forest/wet_drivers/README.md\` on main; NPMS falsification on PR #86.
