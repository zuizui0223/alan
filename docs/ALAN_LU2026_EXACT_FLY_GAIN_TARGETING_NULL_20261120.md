# ALAN experiment: does fly increase preferentially overlap flowers losing nocturnal moth visits?

**2026-11-20 source-verified reanalysis; Draft PR #91**. This is a newly computed conditional-null comparison within Lu et al. (2026) five paired illumination blocks, not a cross-plant species trait effect, pollen delivery, or plant reproduction.

## Original empirical source and estimand

[Lu et al. 2026, Ecology and Evolution, DOI 10.1002/ece3.74368](https://doi.org/10.1002/ece3.74368) published the original Data S1 Excel workbook (source SHA256 \`ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f\`) and R code. The real experiment contained **five independent spatial lit/dark light pairs**. The original spreadsheet has 1,606 plant–pollinator network rows; 1,605 have complete source treatment and diel coding. We discarded the one unlabelled row, without imputing zeros.

The paired source comparison requires exactly the **same plant in the same spatial block with both lit and unlit observations during both day and night**. This yields 33 plant×block cells, from 14 original plant labels in 5 independent lighting pairs. Of these 33, nocturnal Lepidoptera (mostly moth) visits declined in 21.

The original dataset identifies plant labels, number of night moth and day/night fly interactions, and treatment. It does NOT measure same-plant floral pollen received or seeds at these 33 source cells. The 33 cells are not independent randomized treatments.

### Newly tested conditional-null hypothesis

Given the five spatial blocks, is an increase of night or day Diptera *preferentially located on the plants losing nocturnal moth visitation*, more often than expected from the original block-level frequency of each phenomenon?

For each original block b, preserve: n_b = number of paired plant cells; m_b = plants with nocturnal moth visit decline; k_b = plants with fly visit increases in a given diel category. Under the null of fly-gain plant labels exchangeable **within the same block**, overlap X_b is Hypergeometric(n_b,m_b,k_b), E[X_b]=m_b*k_b/n_b. Convolve the five exact hypergeometric probability mass functions to get P(sum_b X_b≥observed) and a probability-ordered two-sided test. This is **not** a randomization of original lamp treatment. It cannot test causal flower selection by flies and does not correct for flower abundance, observer effort or plant identity. It only asks if cross-plant co-occurrence exceeds fixed-marginal expectation.

## Executed source-authentic Results

[Original Excel and exact null GitHub Actions run #44656082048 (attempt 2, SUCCESS)](https://github.com/zuizui0223/alan/actions/runs/44656082048). [Saved reproducible compact output](../results/all_species_forest/biotic_interactions/lu2026_exact_cross_diel_fly_moth_targeting_null_20261120.json).

| Fly gain, same plant and original block | Observed overlap with nocturnal moth loss | Exact block-conditional null mean | Null 95% integer range | Two-sided exact p | One-sided enrichment p |
|---|---:|---:|---:|---:|---:|
| Nighttime | 11 | 11.476 | [9,14] | 1.000 | 0.768 |
| Daytime | 12 | 13.381 | [11,16] | 0.383 | 0.952 |
| Day **or** night | **15** | **16.048** | **[14,18]** | **0.402** | **0.929** |
| Day **and** night | 8 | 8.810 | [6,11] | 0.684 | 0.860 |

Thus, the previously noted **15/21** same-plant fly-gain co-occurrence is NOT enrichment above source-preserved block margins: **15 observed vs 16.05 expected**, with probability-ordering two-sided p=0.402. Original ALAN treatment reallocation of pollinator guild visits is a documented experimental phenomenon, but there is **no support for preferential fly-gain targeting of plants losing moth visits under this conditional ecological test**.

Failure to detect overlap enrichment does not prove absence of temporal pollination compensation: the dataset lacks direct q, compatible stigmatic pollen, seed success, and equal diel opportunity denominators. No cross-time raw counts are summed. Another source paper on a separate focal reproductive experiment cannot be treated as identical block×plant trial without linkage verification.

## Relation to the wet-species observational macro result

This is a real experimental interaction-stage Result, yet does NOT measure whether wet plants are especially vulnerable. Lu's 11 focal species only yield one strict European EIV MOIST match (Geum aleppicum=6.0); the original Tibetan Plateau flora cannot be assigned European ecological indicator values without evidence. There is no same-plant original Lu×Page q linkage, and no formal ALAN→pollinator→wet plant mediation coefficient.

Main document: [Page three static-replacement regimes and null](ALAN_PAGE33_FUNCTIONAL_COMPENSATION_THREE_REGIMES_20261119.md). Keep the wetness macro discovery in the original manuscript and treat this experiment as a candidate mechanism/result about **lack of directed numerical replacement**, not proof of the observed wetness bias.

## Reproduction
- [Original 5-block source-key exact conditional analysis code](../scripts/analysis/alan_lu2026_sameflower_conditional_overlap_20261120.py)
- [Independent hypergeometric mathematics tests](../tests/test_alan_lu2026_sameflower_conditional_overlap_20261120.py)
- [GitHub workflow](../.github/workflows/alan-lu2026-conditional-overlap-20261120.yml)
- [Exact outputs and provenance](../results/all_species_forest/biotic_interactions/lu2026_exact_cross_diel_fly_moth_targeting_null_20261120.json)
