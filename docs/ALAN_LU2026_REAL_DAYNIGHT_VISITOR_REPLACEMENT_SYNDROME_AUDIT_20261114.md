# ALAN — original light experiment: moth-dominant visitor guilds and incomplete fly substitution

**Real original-source analysis SUCCESS (2026-11-14 JST):** [GitHub Actions #43567365631](https://github.com/zuizui0223/alan/actions/runs/43567365631), [aggregate-only artifact #14102987254](https://github.com/zuizui0223/alan/actions/runs/43567365631/artifacts/14102987254). Original primary paper: Lu et al. 2026, Ecology and Evolution, [DOI 10.1002/ece3.74368](https://doi.org/10.1002/ece3.74368), PMCID PMC13645961. This study's ecological phenomenon was already published by its authors; the additions here are a new source-verified syndrome grouping, independently reproducible block-level audit and explicit functional limits.

## Original public data recovered and SHA-pinned

The publisher's genuine supplementary original file is available from [Europe PMC supplemental archives](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles), nested in ECE3-16-e74368-s001.zip. The Excel workbook is actually named Date.xlsx (the original journal Data Availability Statement calls it Data.xlsx); SHA256 **ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f** verified in source runs. The Europe PMC outer ZIP packaging SHA may differ between requests; always reverify the nested original workbook SHA, NOT the dynamic wrapper byte digest.

The author workbook's:
- **visitation** worksheet contains **400** true survey records across 4 survey rounds, 5 experimental blocks, 5 nominal repeat factor levels, illuminated/unlit treatments and day/night. Each diel light arm has 100 source survey records. This is **one** original experiment and **5 field experimental blocks**, not 400 independent regions.
- **species visitation** worksheet contains **1,038** focal plant–visitor rows from **11 actual plant taxa**, both light treatments and day/night. These are a focal subset, not all **29** source night-blooming plant species in the original paper's full network. Plant-by-time absence of a visit is NOT plant absence or independently measured flower availability.
- **ReadMe** specifies five visitor orders and source units are *visits*.

### Source-accuracy flags

The 400 visit records have **3 source duplicate full-factor IDs in nocturnal samples**; some original block×round×treatment cells contain **3–7**, not 5, observation rows. We therefore pair source means at **block×round** (20 matched units per diel), then aggregate study-block effects across five independent blocks. We **do not** fabricate 100 perfectly matched replicate pairs.

The source total abundance is not always exactly the sum of five order columns: **197/200 nocturnal** and **180/200 diurnal** survey records conserve numerically. The other **23 source records** differ by up to 7 recorded visits in one row. These original inconsistencies are reported, not assigned to a chosen insect order. In contrast, all **1,038 species-visitation** source rows individually conserve the sum of the five guilds.

## True source illuminated vs unlit results

Each diel arm has 100 original survey rows per treatment:

| Genuine survey metric | Unlit | Lit | Original lit−unlit | Source proportional change |
| --- | ---: | ---: | ---: | ---: |
| Nocturnal total visits | **1,312** | **848** | −464 | **−35.37%** |
| Nocturnal Lepidoptera visits | **838** | **357** | −481 | **−57.40%** |
| Nocturnal Diptera visits | **54** | **151** | +97 | **+179.63%** |
| Nocturnal Coleoptera visits | 154 | 115 | −39 | −25.32% |
| Diurnal total visits | **1,723** | **1,658** | −65 | **−3.77%** |
| Diurnal Hymenoptera visits | 730 | 584 | −146 | −20.00% |
| Diurnal Diptera visits | **527** | **771** | +244 | **+46.30%** |
| Diurnal Lepidoptera visits | 83 | 40 | −43 | −51.81% |
| Diurnal Coleoptera visits | 167 | 64 | −103 | −61.68% |

These direct source counts independently reproduce the paper's published broad changes; reporting the paper's raw percentages alone is not novel. **The nocturnal Diptera increase of 97 visits is only 20.17% of the 481 lost nocturnal Lepidoptera visit counts.** That is a *visit event count ratio*, NOT 20.17% rescue of actual pollen deposition, seed set, or fruit set, nor a mathematically complete mass balance (source discrepancies noted above). Different insect orders may be dramatically different in effective pollen delivery.

All five original field blocks had negative nocturnal Lepidoptera lit-minus-dark mean differences and all five positive nocturnal Diptera differences. A 2^5 sign-flip *diagnostic*, not verified randomized assignment inference, has a minimum attainable two-sided fraction **0.0625** for each all-same-direction effect. Do NOT promote this to a new independent worldwide p-value: the original paper's own mixed model already established strong within-site visitation changes.

## New source-verified *realized control visitor syndrome* comparison, NOT biological moth dependence

The 11 focal plant species have positive night visitation records under BOTH treatments. We classify a plant as **moth-like by realized visitors** if Lepidoptera represents **more than 50% of all source-recorded nocturnal visits under original unlit control**. This is an explicitly exploratory **visitor-composition grouping** independent of flower color/tube, not a genuine evolved morphologic moth syndrome or effective moth-only pollination.

| 11 focal plant night-use groups | Plant taxa | Original Unlit night visits | Original Lit night visits | Source visit-count change |
| --- | ---: | ---: | ---: | ---: |
| **Lepidoptera-majority night visitors in dark controls** | **5** | **938** | **517** | **−44.88%** |
| Other six documented night visitors | **6** | **249** | **231** | **−7.23%** |

The five majority-class plants lost **425** documented nocturnal Lepidoptera visits (680→255), while gaining 46 Diptera visits (87→133). The six non-majority plants lost 22 nocturnal Lepidoptera visits (95→73) but gained 29 Diptera visits (19→48). Thus **5/5** originally Lepidoptera-majority focal plants showed fewer lit-site nocturnal Lepidoptera visits; **6/11** focal plants showed both moth visits decreasing and fly visits increasing.

This suggests the **identity and degree of dependence on light-sensitive visitor channels**, rather than a simple binary moth-visited label, can help predict the magnitude of source-observed visitation changes. However the 5-versus-6 difference is **within one original experimental research system, not a replicated causal cross-syndrome meta-analysis**. Control-majority group assignment and treatment differences share one original data source; flower availability, baseline frequency and uneven species detection confound the comparison. No p-value is reported for this post-hoc 5-vs-6 grouping, and the fact that all six group comparisons are positive-record restricted prevents treating undocumented local guilds as absent.

## Macroecological extension and limitations

What is newly useful for global ALAN research is a reproducible, source-backed operational definition of **realized nocturnal visitor-channel reliance and alternative-visitor response**, rather than treating a globally recorded flower shape or any GloBI non-record as pollinator dependence. This must then be validated against true pollen deposition and seed set, locally measured flowering-time overlap, selfing/pollen limitation and predation.

Critically these data do **NOT** contain actual seed/fruit outcomes, stigma pollen, predator attacks, predator substitute pollination, or experimental island/mainland field contrasts. The Tibetan meadow is one terrestrial site. The Swiss ALAN phenology Mendeley dataset [10.17632/xjmn5jhxyn.2](https://doi.org/10.17632/xjmn5jhxyn.2) could not yet be retrieved via its public metadata API (3 candidate endpoints returned HTTP 401); no Swiss mismatch coefficient has been generated. An earlier 562-island photo–ALAN evidence registry with 52 actual independently documented night/day reproductive species and incomplete syndrome trait coverage remains **a separate observational source**.

## Reproduce original source

Original outer/inner supplementary ZIP is read at GitHub Actions runtime and the nested XLSX content is SHA-verified. Python analysis: scripts/analysis/alan_original_Lu2026_daynight_guild_syndrome_source_20261114.py. Five-block and source row-conservation tests: tests/test_alan_original_Lu2026_daynight_guild_syndrome_source_20261114.py. Real workflow: .github/workflows/alan-Lu2026-original-5block-syndrome-visitation-20261114.yml. All resulting uploaded files are **aggregate-only**; no original plant×block×repeat row-level workbook or copyright source is published in this fork. Prior results/PR #64 and main remain unchanged.
