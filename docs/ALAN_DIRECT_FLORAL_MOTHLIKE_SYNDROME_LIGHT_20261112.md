# ALAN × species-direct moth-like floral phenotype: real source test (2026-11-12)

**Status:** Original 2019/21/23/24/25 plant-photo + same-year nightlight replay [#43161009513](https://github.com/zuizui0223/alan/actions/runs/43161009513) SUCCESS; a stricter fail-closed geographic-sample check is included in this branch. **All results exploratory.** This is a morphological resemblance study, **NOT verified nocturnal/effective moth pollination**.

## Question and fixed empirical source

Are GloBI old-only 713 species with a negative source-photography light component enriched for independently recorded white, long-tubed, deep-throated flowers?

- ALAN original **10,698** 0.05° photo-cell×year records (**933 island, 9,765 mainland**), WCVP accepted species **8,581**, and **78** 5° geographic clusters. Same photo counts, weather, population, 2019/2021/2023/2024/2025 four sampled monthly VIIRS values, year×region FE and island/mainland balanced weights. Recomputed GloBI baseline β=**−1.901720622**.
- Independent upstream observed species-DIRECT source: zuizui0223/island run **34191508045**, artifact **10042377593**, direct_species_trait_ledger.csv.gz SHA-256 **7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577**. High/medium recorded trait states, no genus or LLM proxy. Only unambiguous single-state **white** and one of **tubular / salverform / funnel_trumpet** counts as BROAD moth-like; STRICT additionally requires source-direct **deep**. Separate high-only sensitivity.
- Exact accepted scientific species-name joins only; multistate, unrecorded, and unresolved traits remain **UNKNOWN**, never "not moth-like". Source restrictions prohibit redistributing original species × trait or species × photo-site identities; this branch exports aggregate-only models.

## Actual ALAN photographed species-level morphology coverage

| GloBI label group | All accepted photographed species | Both single-state color and form known | White + tube-like broad | White + tube-like + deep |
|---|---:|---:|---:|---:|
| Old-only | **713** | **128** | **5** | **0** |
| Recorded by both | **612** | **91** | **7** | **0** |
| New-only | **132** | **27** | **1** | **0** |
| Neither recorded | **7,124** | **684** | **88** | **7** |
| **All** | **8,581** | **930** | **101** | **7** |

Among taxa with both measured color and form, only **5/128=3.9%** of the old-only 713 are broad moth-like, against **7/91=7.7%** among positives retained in both GloBI versions. Fisher exact odds ratio **0.488**, nominal two-sided **p=.243**, not an enrichment. ALL 7 strict white+tubular+deep photographed plant species belong to "neither GloBI label recorded"; this is NOT proof of true moth absence. High-only sensitivity: **4/84** old-only vs **1/49** common, Fisher p=.651; no supported result.

## Frozen site×year source ALAN regressions (exploratory)

Response is the percentage of distinct observed accepted species per source photo-cell×year with an independently documented species-direct morphology attribute. Each β is a spatial **island minus mainland × log(1+same-year VIIRS)** contrast, in percentage points per log brightness, not a local pollen-deposition or causal effect.

| Photo-recorded endpoint | Positive island / mainland photo-cell×years | Independent 5° geography blocks with both recorded positive sites | β | nominal 95% cluster CI | p |
|---|---:|---:|---:|---|---:|
| Broad white + tube-like morphology (101 source species) | **52 / 518** | **15** | **+0.2891** | **[−0.1298,+0.7079]** | .173 |
| Single-state white flower | 450 / 4,462 | 56 | +0.5063 | [−1.9016,+2.9141] | .677 |
| Tube-like flower form (regardless of color) | 338 / 3,701 | 56 | −0.2818 | [−1.6423,+1.0786] | .681 |
| Strict white+tube+deep (7 species) | **5 / 41** | **2** | **NOT ESTIMABLE WITH SUFFICIENT INDEPENDENT GEOGRAPHY** | — | — |

An independent post-hoc matched response among *only known color+form photo species* (at least two classified species per eligible cell-year, original same-year geographical matching) has **β=−0.9114 [−4.2073,+2.3845], p=.572** (32 positive island and 237 positive mainland cell-years across 11 common-positive blocks). No robust negative morphology-specific association.

**Source High-only sensitivity** broad morphology across all original photographed species: **β=−0.0012 [−0.2550,+0.2526], p=.992**, only 30 positive island and 199 positive mainland cell-years with 10 shared-positive geographic blocks. Neither effect is statistically resolved.

## Support guardrails and scientific claims

GloBI **old-only × white+tube** has only **5** positive island photo-cell×years and **2** shared 5° blocks; GloBI common × white+tube has 5 positive island cell-years and 3 shared blocks. The final analysis refuses to report pseudo-precise subgroup effects where fewer than 20 positive island/mainland photo-cell×years or fewer than 10 shared-positive geographical blocks. An unbalanced spatial result or modeled coefficient with inadequate cross-land support is not a demonstrated phenotype effect.

**Interpretation:** independent species-direct floral evidence DOES NOT support a shared, distinctive white/tubular/deep moth-like morphology in 713 old-only plants, nor a universal negative light response of morphologically moth-like plants. This is **not proof of ecological absence**, because only **930/8,581** observed source species have usable paired color+form, and direct depth coverage is smaller. Recorded morphology is not observed night pollination; all flower color, shape, deep tube and group calls are species-level external traits, not sample-site floral measurements.

Next: independent actual nocturnal pollen deposition / day-night exclusions, wild-native island establishment checks, stronger source-approved species-direct floral morphology, and genuinely repeated plant reproductive success across nighttime-light changes. Never use absent GloBI relationships as negative biological controls or convert photos of ornamental plants into native-flora prevalence.

## Sources, audit and reproducibility

- [Real source Actions #43161009513](https://github.com/zuizui0223/alan/actions/runs/43161009513)
- [Final stronger evidence-support workflow](https://github.com/zuizui0223/alan/actions/workflows/alan-direct-floral-syndrome-ALAN.yml)
- Source key SHA: 37f3c2c519871b6905ffa81c5a14536b248611861ad143b1d4d77d95b7c6b00d
- Scripts: scripts/analysis/alan_direct_flower_syndrome_light_20261112.py; tests/test_alan_direct_flower_syndrome_light_20261112.py
- This Draft PR is stacked on the frozen old-v0.11 attribution PR #43; main and previous results unchanged.
