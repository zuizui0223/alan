# ALAN B: Silene latifolia source-light contamination in nocturnal functional dependence (2026-11-15)

## Source and actual executed result

The existing file results/night_exclusion/experiments.csv contains **42 original day-only/open comparisons and alternate endpoint rows** across pollination-exclusion publications. results/night_exclusion/species.csv stores a study-equal mean night-access loss for 12 plants. These are day/night access data, **NOT the same as plant fruit/seed response to experimental ALAN**.

For the A hypothesis a source contrast must NOT already blend experimental ALAN assignments inside its proposed dark daytime/open V. **Macgregor et al. 2019, Ecosphere DOI 10.1002/ecs2.2550** appears in the old night-exclusion table with D=0.3935018051, O=0.7346570397, so V=0.4643734644. Its Figure 2 caging estimates were **pooled across the five distinct lighting treatments**, not obtained exclusively from the unlit plots. The original Results also report pooled approximate day-only 37.6%, open 72.0%, night-only 52.8% ([Wiley original article](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.2550)). Therefore these numbers cannot serve as an unexposed baseline moderator for the same experiment's own ALAN outcome.

## Verified sensitivity, historical results preserved

[Actions #43706444027 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43706444027): six original-source audit tests passed. Artifact #14158119852, SHA256 of ZIP 3aff1a454d0196dcf6b822dd836af3859fb44caa88229b9691e4d3bc56084dbb. Script scripts/analysis/alan_Silene_lat_dark_only_V_gate_20261115.py reproduces the original 12 plant species means from source; then creates a **new separate strict sensitivity** excluding Macgregor's mixed-light row. No original data or previous published results overwritten.

| Silene latifolia source summary | Independent source studies | Study-equal V = 1-D/O |
|---|---:|---:|
| Original source summary, with mixed-light Macgregor | 5 | **0.6718269902** |
| Strict source sensitivity, omitting pooled-light Macgregor | 4 | **0.7236903716** |
| Strict minus original | one study removed | **+0.0518633814** |

The remaining **6 used source contrasts** from 4 papers (some different site/year/metric arms) span **0.5057–1.0000**. Individual study means: Anderson_2015 0.9146 (two Saskatchewan years); Bartelmess_2006 0.8753 (two Tennessee distances); Scopece_2018 0.5992 (Italy); Young_2002 0.5057 (Colorado). Heterogeneous locality, year, flowers and reproductive denominator make a fixed species-wide V an uncertain proxy; this is NOT within-site effect variation, adaptation or evolution demonstrated by the new audit.

Other eleven species means in the old appendix are unchanged by this specifically targeted correction. Ambient skyglow at the four remaining non-ALAN source sites is not measured uniformly; the label is **"no intentionally assigned ALAN treatment included"** rather than physically verified 0-lux night.

## Why this matters for A and B

Nighttime functional dependence is not a species-immutable "moth pollination syndrome" trait. For the same UK population as Macgregor's ALAN experiment, we require **dark-only** diurnal-only and unrestricted pollinator-exclusion outcomes at original lamp-site and round, together with clustered uncertainty. Borrowing another country's species mean (even corrected to 0.724) is external transport with unmodelled environmental and demographic variation, not independent local moderator data.

The 2025 small-fruit 850-flower Figshare study was **already represented by Rubus idaeus** in the original 42-row source file; the newly extracted five species×year contrasts are an expansion of the very same publication, NOT a new second source study. The island×ALAN photograph model was not refitted; its response is not experimental fruit or seed set.

## Reproduction and limits

Audit outputs:
- results/night_exclusion/S_latifolia_dark_only_source_gate_20261115.json
- results/night_exclusion/strict_night_V_by_species_20261115.csv
- tests/test_alan_Silene_lat_dark_only_V_gate_20261115.py

Unchanged historical inputs: results/night_exclusion/experiments.csv; results/night_exclusion/species.csv.

A genuinely independent local dark-only V (with source variance/covariance), an open-lit/open-dark seed/fruit rate and its true treatment installation cluster remain unavailable for meta-regression. Confirmatory beta_V and island×ALAN reproductive interaction remain **unestimated**.


## 12-species observational light-gap sensitivity (post hoc diagnostic only)

This pipeline additionally re-evaluates the **simple species-level ordinary least-squares descriptive slope** of the old 12 species' reported island-minus-mainland photographic night-light association gap (the old species.csv gap column) versus the source-derived nighttime dependence V. It is a **different dependent variable from plant reproduction** and is not an experimental ALAN response.

| All 12 species, gap OLS coefficient per +10 percentage points in V | Estimated point slope |
|---|---:|
| Original V (including the mixed-light Macgregor source) | **+0.002326** |
| Strict source-based V (remove Macgregor source from S. latifolia) | **+0.002230** |

The point slope changes by roughly **−4.1%**, retaining a positive sign, rather than the proposed negative relationship. This is NOT a new source-backed robust significance verdict: the original 78-region clustered uncertainty, spatial sampling, species-composition changes and original photographic data were not refitted and the original reported p-value **may not be reused as a p-value for the corrected slope**. It is only a sensitivity of the 12-species coefficient conditional on the original frozen species gap values.

Again: an observational species-composition photographic coefficient cannot be called an ALAN-caused change in seeds/fruit or an ALAN×nighttime dependence plant reproductive meta-effect.
