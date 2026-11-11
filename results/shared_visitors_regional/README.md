# Regional flower visitors in the five-year island–mainland comparison

The frozen estimate is −1.901721 percentage points per unit log1p(VIIRS) for the island minus mainland slope in the share of photographed plant species with a documented moth-family flower-visit link. The frozen sample has 10,698 site-years across 78 five-degree cells. The regional analyses below describe a selected subset.

## Coverage

We screened [UK PoMS FIT Counts](https://doi.org/10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618), [EuPPollNet](https://github.com/JoseBSL/EuPPollNet/tree/31716d21d5aad8c27c1d91dd0620149ef42e8836), [CropPol](https://github.com/ibartomeus/CropPol/tree/99c69dbab4fa7cbdb630ff3e9bd9989b3de9548b), [PollAgri](https://doi.org/10.48579/PRO/YV7MXG) and [HJ Andrews meadow surveys](https://andrewsforest.oregonstate.edu/data/datacatalog/SA026). The latter three contributed no exact plant × cell × year matches to the frozen frame. Working taxonomic labels required a binomial matching a frozen accepted name or a WCVP-resolved name. Aggregate or unresolved labels provide no species-level support.

Of 544 plants observed on both land classes somewhere in the frozen frame, 38 had positive survey evidence on the same species in the same cell and calendar year as a photograph. The 42 plant–cell–year combinations span six cell-years (2019, 2021 and 2023). There were no matches for 2024 or 2025. Positive species counts by survey region were 7 (eastern Britain), 31 (Germany), 2 (southern Sweden) and 1 (Menorca), with overlap between species lists. Only one combination had photographs of the same plant on both land classes in the matched cell-year. No matched cell-year had positive visitor survey records from both land classes according to the frozen GSHHG classification. Missing coverage cannot be interpreted as zero visits.

## Regional models

The response for the shared-plant analysis is the percentage of all photographed plant species per site-year belonging to the 544. A parallel fit uses the original moth-linked plant share. The model includes region-by-year fixed effects, the frozen covariates, within-stratum land-class balancing, and interactions between log1p light, island class and the regional visitor metric. A region-wide visitor metric is assigned to both land classes. Region-years with at least two light values on each side and strictly overlapping light ranges define the reported comparison. All site-years in such strata enter the slope estimates.

The UK absolute-frequency metric is the mean count of bees/wasps, flies and beetles per 10-minute focal-flower count, including sampled zeros (49 counts, 490 minutes across all three years). Plant means receive equal weight within a cell-year. Other protocols lack a comparable per-plant visit-rate denominator. Across studies, the composition metric is the expected number of these three orders represented in two randomly drawn visitor individuals from the observed pool, including other visitors in the denominator. Source-by-plant estimates require at least two individuals. This measures order representation and cannot estimate total visits or pollen delivery.

| Regional predictor | Cells / cell-years with light overlap | Island / mainland site-years | Shared-544 slope gap at predictor Q25 → Q75 | Original response gap at Q25 → Q75 |
|---|---:|---:|---:|---:|
| Visitors per 10 minutes | 1 / 2 | 22 / 18 | +15.82 → +11.59 | +21.50 → +12.98 |
| Representation of three visitor orders | 3 / 4 | 31 / 72 | +6.96 → +9.75 | +6.48 → +11.56 |

Both response columns use percentage points per log1p light. Q25 and Q75 are computed over unique eligible cell-years. The absolute-rate comparison uses between-year variation within a single region. The cross-region composition comparison pools unequal sets of plants and protocols. Its estimated gap widens as order representation increases. Neither analysis identifies a global attenuation of the frozen negative slope difference or a compensatory pollination mechanism. No cluster-based confidence interval was reported because fewer than 10 independent regions qualified.

## Reproduction

Run scripts/analysis/alan_shared_visitors_regional.py with positional arguments for the existing v0.11 environmental/WCVP inputs (including coast), MOIS signed directory, downloaded survey directory (filenames and hashes in source_sha256.csv), and output directory. The frozen target digest and coefficient are checked before joining. regional_plant_support.csv lists positive matches, shared_544.csv accounts for every candidate, survey_land_support.csv records survey-side assignments, the regional CSVs contain the actual predictor values and light eligibility, and summary.json contains all model estimates.
