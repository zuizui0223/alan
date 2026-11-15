# ALAN B | Global bloom-time × tube-shape syndrome versus functional nocturnal dependence (2026-11-15)

## Authentic source and design

Kendall & Nicholson (2025) Ecology Letters DOI 10.1111/ele.70036. Both genuine author source files were fetched from the fixed GitHub commit 1945b7c114ad87000608ce6018d0280001c04ef1 and SHA256 validated:

- Original original DvN reproduction and exclusion file: 1,652 treatments, SHA256 1d9c4c778206584130308b99ba21d902b0ccf1047118dfcf52abaec9d652e2aa.
- Original floral trait file: 142 taxa, SHA256 9c4a00d74c954eb448330066071ee0c584a9df66d96911c011cef5cd1afc7173.

Within-source study, taxon, location, year, time, treatment, reproductive endpoint and unit matching produces **340 D/O reproductive contrasts from 113 original study IDs**; these are **not** 340 artificial light (ALAN) experimental reproduction observations. Functional nighttime-access marginal reliance is V = 1 - D/O (dark day-only divided by open), but ambient light is not known to be physically 0 lux. The author trait table has 69 night-blooming species, 30 day-blooming, 19 both, 24 unknown, of which 87 are tube shaped. Trait matches among the full 340 reproductive contrasts are night 210, day 60, both 46, unknown 24; tube 253, non-tube 85, unknown 2. All unknown values remain in the complete source set.

**Critical source correction:** the Macgregor_2019 single floral reproductive comparison uses day/open means pooled ACROSS its 5 experimental lighting conditions. This is invalid as a light-free baseline moderator in A. Both source 340 data and a STRICT sensitivity with this light-pooled source removed (339 comparisons, 112 studies) are kept; no original rows deleted.

## The original study-equal descriptive comparisons

Within each paper, source contrasts are represented by their median; these paper medians are then equally weighted and compared. We do not mix fruit formation with seed per successfully formed fruit or seed per ovule. The contrast is night blooming + tubular versus day blooming + non-tubular flowers.

| Reproductive denominator | Night-tube V: source comparisons / source studies | Day-nontube V: comparisons / studies | Descriptive difference | Publication bootstrap 2.5–97.5% interval |
|---|---|---|---:|---|
| Fruiting probability, all 340 source records | **0.5336** (51 / 22) | **0.3021** (8 / 6) | **+0.2315** | −0.4598 to +0.8103 |
| Fruiting probability, STRICT source with Macgregor removed | **0.6029** (50 / 21) | **0.3021** (8 / 6) | **+0.3007** | −0.4715 to +0.8576 |
| Seed count per formed fruit | **0.2740** (66 / 24) | **0.1434** (7 / 4) | **+0.1307** | −0.2282 to +0.4931 |
| Seed per ovule | **0.4271** (2 / 2) | **0.5040** (2 / 2) | −0.0769 | Suppressed, too few papers |

The meaning of 0.60 is: the functional reproductive metric in the day-only access treatment was approximately 60% lower than in the day+night unrestricted access treatment **as a within-source study-equal median**, NOT an expected 60% seed loss under ALAN. A 0.30 difference here means 30 percentage points in the relative-loss index, not 30% absolute fruit set or a comparative light impact.

**There are ZERO papers in these primary subtype comparisons with both the night-tube and day-non-tube source groups in the same original study.** The apparent positive syndrome pattern is entirely across-study/taxon/geography, and cannot be turned into a causal light sensitivity or a verified nocturnal moth pollen service effect. Day-blooming tubular plants have only ONE original study in the fruiting-probability layer and night-blooming non-tubes only TWO; the effects of blooming time versus shape cannot be separately identified as a reliable factorial interaction.

### Whole-publication leave-one-out source robustness

- STRICT fruiting probability: delete each of the 21 night-tube or six day-nontube source papers in turn (27 removals); all 27 differences stay positive, range **+0.0377 to +0.5638**. The independent paper-bootstrap distribution nevertheless crosses zero.
- Original fruiting probability including Macgregor mixed light: 25/28 deletions positive, three negative, range **−0.0316 to +0.4946**.
- Seed-per-formed-fruit: 28/28 deletion differences positive, range **+0.0026 to +0.2587**, but original paper bootstrap crosses zero.

Twelve thousand independent draws of publication IDs per original group give descriptive intervals above. They are **not causal or phylogeny/within-plant variance-corrected confidence intervals**. Paper-level resampling cannot correct publication selection or same research group / related taxa across different papers. The much smaller daytime comparison group is a particularly strong limitation.

## Ecological inference (restricted and falsifiable)

The exploratory direction is compatible with night-blooming tubular plants having greater realized nighttime-access benefit in reproducing. But **blooming time and tube morphology are crude traits, not verification of nocturnal insect effective pollination or their specific ALAN susceptibility**. The inequality can also arise from study geography, sampling, pollinator assemblages, and outcome denominators.

Use this B result as a descriptive guide to stratify additional primary studies. Do NOT generate a universal species-level V solely from a morphology tag, or treat these 340 source comparisons as experimental ALAN fertility effects. The A confirmatory meta-regression beta_V and the predicted island×ALAN fertility interaction remain **UNESTIMATED** because aligned reproductive lighting effects, local independently measured V and source light-cluster covariance are not recoverable.

## Reproducibility

GitHub Actions success logs re-fetch and hash-check both primary author files; seven source-and-leave-one-paper unit tests protect endpoints, actual study IDs, the Macgregor exclusion and no fabricated experimental causal inference.

- Code: scripts/analysis/alan_B_P1_original_syndrome_functional_20261115.py
- Tests: tests/test_alan_B_P1_original_syndrome_functional_20261115.py
- Workflow: .github/workflows/alan-B-P1-340-global-floral-syndrome-20261115.yml
- Result artifact: results/alan_B_P1_global_floral_syndrome_functional_v_20261115.json

Original 1652 source rows and 142 taxon trait rows are NOT redistributed in this project output; only derived aggregate values/metadata are retained. Historic previous results are not overwritten.
