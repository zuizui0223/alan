# ALAN湿潤種偏りの機構監査：真の草地内・調査回内の花被度比較と観測選択の影響

2026-11-20 | original Giavi et al. 2021 field source | **SOURCE-SUPPORTED DESCRIPTIVE ONLY** | main unchanged

## Why this falsification is necessary

Our large observational macro result (+0.534 adjusted indicator-unit moisture difference of illuminated-site recording-decline candidates) is retained as a *finding about plant trait representation*, not a causal ALAN effect.

Earlier independent-source Giavi 6-meadow ALAN data show a nominal negative moisture–floral-cover difference after same-sampling-event FE: PLANTATT F beta −0.769 (descriptive six-site whole-meadow label signflip p=.03125), Döring F beta −0.646 (p=.25). Both survive removing any SINGLE plant species, but **only two truly separate physical meadow sites jointly contain a wet F>7 and a drier F<=7 named plant**. The earlier 53/77 matched plant–event rows cannot be used as 53/77 independently illuminated meadows.

The original authors measure **area covered by flowers** along their 100 m transects. Their flower surface proxy (plab) is measured using discrete approximately 5 cm diameter circles shortly before visitor sampling. This is NOT a flower count, a species' opening/closing time, a seed outcome, or a pretreatment baseline. Original paper [Giavi et al. 2021, Nat Commun](https://doi.org/10.1038/s41467-021-22011-8), authentic original [Zenodo Data RData](https://doi.org/10.5281/zenodo.4540407) MD5 f60b21d82c163cbb53520fda0bc61985.

## New original-source within-meadow-and-exact-sampling-event comparison

[Successful original source code and synthetic tests, CI #44670559432](https://github.com/zuizui0223/alan/actions/runs/44670559432). [Source exact JSON and per-site results](../results/all_species_forest/biotic_interactions/giavi2021_original_meadow_event_wetdry_flower_coverage_20261120.json).

For each SAME physical original meadow AND SAME one of 37 true original paired sampling events, take the mean log(1+flower surface proxy) for the observed wet plant species (F>7) minus the mean for observed drier species (F<=7), separately for lit and dark. Difference those two same-event group gaps: [(wet−dry)lit − (wet−dry)dark]. A source plant must exist in both actual original arms, and every species within an event receives equal weight. Every physical site's eligible event gaps are averaged, then the source-level estimate is the mean of the independent meadow means.

| Empirical support (all from original RData) | UK PLANTATT F | UK Döring F |
|---|---:|---:|
| Exact same original meadow×sampling event with both wet/drier plants | **3 events** | **7 events** |
| Independent original physical meadow pairs supporting comparison | **2** | **2** |
| Negative per-meadow paired wet-drier floral gap | **2/2** | **2/2** |
| Mean of two original meadow gaps (log1p floral-cover unit) | −3.346 | −0.895 |
| Original meadow 1 wet−dry cover Δ | −0.832 | −1.297 |
| Original meadow 6 wet−dry cover Δ | **−5.860** (only 1 eligible event) | −0.493 |
| Smallest two-sided whole-meadow signflip p, even if directions perfectly concordant | **0.500** | **0.500** |

The exceptionally negative one-event PLANTATT value at original meadow #6 is NOT a reliable species-wide estimate. At the meadow with 3+ original qualifying events, Döring estimates are far less extreme; comparing these approaches as independent experiments would be wrong.

**Interpretation**: two true meadow contrasts lean negative, but the actual design is incapable of independent population-level wetness×illumination inference, regardless of apparent statistical significance computed from 53/77 plant-event rows.

### The decisive selection check: which plant measurements have TWO original arms?

The authors' source includes a plant×site×sampling-event in just ONE arm for many cases, and a true paired lit/dark measurement for others. An unreported arm is **not zero floral cover**, extinction, or failure to flower.

Döring F original data:
| Plant ecological affinity | Both C/T arms recorded in exact source event | Only a single recorded arm | Both-arm source coverage |
|---|---:|---:|---:|
| Wet F>7 (2 distinct plants) | **8 of 27** original plant×event keys | **19 of 27** | **29.6%** |
| Drier F<=7 (6 distinct plants) | **69 of 130** | **61 of 130** | **53.1%** |

PLANTATT F scores (fewer exactly scored plant names): wet still **8/27** both arms, drier **45/102** both arms (44.1%). These are raw repeated plant×event coverage fractions, not independent Bernoulli trials for p-values. The lower paired coverage among wet plants suggests potential **selection into the analysable matched cohort** that was not addressed by the original site-conditional field trait slope.

### Crucial dry-only negative control: the negative moisture trend is NOT restricted to F>7 plants

Restrict entirely to original plant species with F<=7, removing BOTH wet F=8 species, and retain the same sampling event FE for the floral resource contrast:

| Dry-only source subset | Original plant labels | Original paired plant-events | Within-event F slope on lit-minus-dark log floral cover |
|---|---:|---:|---:|
| PLANTATT F 4–5 | **5** | **45** | **−1.009** |
| Döring F 4–7 | **6** | **69** | **−0.637** |

These are descriptive coefficients, not independently inferred biological effects; only small F score ranges remain. They show that the negative continuous-F slope can arise **entirely among plants not classified as wet (F>7)**. Thus, calling this an independently verified special sensitivity of wet plants is not justified by this source.

Site-by-site F slopes can be positive or negative; even among the six original meadows only two contain any scored wet+dry comparison. The original experimental light assignment has not been verified as randomly assigned in our raw source and pretreatment flower resource was not measured. Neither whole-site label flips nor this source selection check establishes an ALAN causal wetness effect.

## Consequence for the manuscript

Keep macro moist-plant bias Results 1 unchanged. For the mechanisms chapter:

1. Source Giavi 2021 empirically measures ALAN lit/dark *flower surface proxy* and daytime insect visits in a 6-pair field setup. The original paper's genuine insect guild redistribution is experimental evidence.
2. Our newly computed measured floral proxy contrast suggests a quantitative **flower availability vs visitor allocation distinction**. But this is NOT a robust wet-only ecological mechanism due to poor within-meadow overlap, strong source arm selection, and a similarly negative within-drier F gradient.
3. Page original 72 pollination networks provide independently measured per-visit efficiency heterogeneity, but no same-species ALAN causal stage and no wetness moderator. Lu 2026 5-block fly increases do not preferentially target plants losing moth visits after exact fixed-marginal block null (15 observed vs 16.05 expected, p=.402).
4. More compelling wet-plant causal mechanism experiments must have >=3 physical independently treated source site pairs each containing >=2 wet and >=2 drier plants on the SAME sampling date, pre-light floral cover or randomization ledger, as well as compatible pollen/seed or growth/fate endpoints; pooled multiple papers from the same light installation count as one experiment.

The best current Results phrasing is:
> In an independent six-paired-meadow field ALAN experiment, floral resource availability varied with the ecological moisture score of flowering plants, but only two original meadow pairs supported direct within-meadow wet-versus-drier comparisons. Focusing on same-meadow sampling events and inspecting original one-arm recording showed that the evidence cannot isolate a wet-specific treatment response.

**Do not promote the new source-level figures as independently confirmed causes of the macro wet-plant bias, and do not call observational floral cover association a causal mediation effect.**

## Reproduction and source provenance

- [Original same meadow, same sampling event and missing-arm source audit script](../scripts/analysis/alan_giavi2021_same_meadow_wetdry_floral_20261120.py)
- [Synthetic exact event arithmetic and true physical-site evidence tests](../tests/test_alan_giavi2021_same_meadow_wetdry_floral_20261120.py)
- [Workflow source replay](../.github/workflows/alan-giavi-same-meadow-wetdry-20261120.yml)
- [Verified original source aggregate](../results/all_species_forest/biotic_interactions/giavi2021_original_meadow_event_wetdry_flower_coverage_20261120.json)
- [Prior original 6-meadow source flower-FWL analysis and sensitivity](ALAN_GIAVI2021_TRUE6MEADOW_FLORAL_AVAILABILITY_WETNESS_FWL_20261120.md)

Original independent Kew BIFloraExplorer UK indicator csv Git blob SHA1 e6a1da929af81f370d2a1b79177b489e3e573a06.
