# Lockett2022 independent design and allocation audit

Reviewed 2026-11-20. Verdict: original observations recovered, but conventional independent-arm extraction remains on hold.

## Primary-source summary (under 200 words)

[Primary article](https://doi.org/10.1016/j.envpol.2022.119803) describes 24 light-screened plots in Kongwak, Australia (−38.513677, 145.708209), assigned alternating ALAN/control in a checkerboard, rather than randomized lighting allocation. Six potted Eucalyptus camaldulensis trees per plot were subsamples. Psyllids were introduced in mesh bags; colonized trees moved after seven days, including within-treatment moves, producing four exposure histories. First treatment covers prior tree conditioning and parental oviposition; second treatment covers subsequent juvenile development. These are not interchangeable instantaneous and mediated effects.

At week 40, 23 trees with constant lighting history were harvested, reportedly one per final plot (12 ALAN, 11 control). Shared initial environments can remain relevant despite separate final plots. Biomass allocation and morphology are distinct from total biomass. Reported increased lerp production does not establish increased consumption, predation, or plant mortality. Failed oviposition samples were excluded from subsequent-stage analyses; those stage-conditional outcomes cannot be read as unconditional survival. The article reports comparatively limited effects on other psyllid life stages.

## Independent original-data audit

Read all six sheets of the tree workbook and metadata, and both CSV schemas in [Mendeley version 1](https://data.mendeley.com/datasets/pfrrkgs878/1). Parent verified download hashes; this audit checked data contents, not a second full hash run. Joins below use unique tree IDs and the week-38 Growth records. No source values changed.

### Harvest allocation and conflicting labels

The biomass sheet contains 23 unique trees, each with a distinct final `plot`; all 23 final plot numbers agree with Growth `plot38`. However, they came from only **18 distinct initial plots**, with five shared-origin pairs:

| Initial plot | Harvested trees |
|---|---|
| 1 | 5, 6 |
| 2 | 9, 7 |
| 3 | 14, 15 |
| 13 | 76, 73 |
| 16 | 91, 96 |

Consequently 23 terminal trees cannot simply be assumed mutually independent over the complete exposure history. All retain the same lighting treatment according to Growth, but changed physical plots. Separate destination plots do not remove dependence induced by a shared earlier plot. Conversely, collapsing solely by origin discards destination dependencies, including a physical plot functioning in both periods. A declared exposure-history model or a justified conservative analysis is needed; do not substitute an effective sample size of 18 without deriving its covariance.

Four biomass-sheet treatment labels disagree with Growth and PCA:

| Tree | Biomass sheet | Growth history | PCA treatment | Psyllid corroboration |
|---|---|---|---|---|
| 91 | Dark | light_light | ALAN | No row |
| 14 | Dark | light_light | ALAN | 14N: Light→Light, plots 3→7 |
| 7 | Light | dark_dark | No ALAN | No row |
| 60 | Light | dark_dark | No ALAN | 60N/S: Dark→Dark, plots 10→8 |

These swaps preserve 12/11 arm totals, so a sample-size check misses the conflict. Multiple-source agreement supports a reconciliation candidate, not silent correction. Biomass East/West labels are reversed relative to Growth for the joined trees. PCA uses a differing Northness convention. Preserve each source's position fields until mapping is documented.

### Psyllid outcomes and denominators

The main CSV has 182 bag records and identifies both initial and final plots. `Use4ALAN=Yes` selects 117 bags:

| First→second exposure | Bags | Trees | Initial plots | Final plots |
|---|---:|---:|---:|---:|
| Dark→Dark | 30 | 19 | 11 | 11 |
| Dark→Light | 30 | 20 | 10 | 12 |
| Light→Dark | 30 | 19 | 12 | 10 |
| Light→Light | 27 | 16 | 11 | 10 |

These are observed allocation counts, not independent lighting n. Two bags can share a tree, trees share initial/final plots, and different histories share plots; four cell means need covariance. All observed initial-plot and final-plot treatment mappings are internally single-valued. The subset contains 104 `FAILED=NO` and 13 `REZERO`; all 32 `FAILED=YES` rows have missing `Use4ALAN`. Hence this flag alone cannot define the denominator for oviposition failure or parental survival. Each endpoint needs its own explicit inclusion rule, timing and at-risk denominator.

Metadata defines `PROP_ESTAB_STD1` as division by the observed maximum and `PROP_ESTAB_STD2` as capping values at one. Neither is an unmodified biological probability. Metadata warns that `LERPPROD` differs from the published analysis numerator/denominator: use `lerp / LERP_ESTAB_FINAL` for that estimand, after checking validity; do not directly treat the archived LERPPROD column as identical.

### Extraction decision

Original harvest biomass and allocation components are numerically present, as are stage-specific psyllid counts and plot histories. Their availability supports further reproducible reconciliation, not immediate independent-arm meta-analysis. Hold pooled effect-size entry until treatment conflicts, crossed exposure covariance, selection and endpoint denominators are resolved. Report author-described morphology and lerp findings as published evidence, clearly separate from self-extracted effects. No attack rate or predator-mediated effect was measured here.
