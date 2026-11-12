# ALAN: Are the 713 plants positive only in the old GloBI cache a coherent ecological guild?

**Date:** 2026-11-12  
**Status:** original signed photo/taxonomy/GloBI/VIIRS/WorldPop source replay SUCCESS. Exploratory *documentation-ecology diagnostics*, not a test of local effective moth pollination.

## Frozen source and exact coefficient decomposition

Original observations and environmental model held fixed:
- 2019, 2021, 2023, 2024 and 2025; 10,698 cell×years (933 island; 9,765 mainland), 78 independent 5° spatial blocks and 195 year×block fixed-effects groups.
- 8,581 photographed WCVP accepted species; each accepted species counted once per source cell×year.
- Old published GloBI moth-positive plant label versus full v0.11 accepted-ID label. No photos, light pixels, denominators or covariates are swapped.
- Actual original scientific replay: [GitHub Actions #43157177667](https://github.com/zuizui0223/alan/actions/runs/43157177667) **SUCCESS**, artifact ID **13916655585**. Test suite and signed input checksum gate passed; new results held in downloadable *aggregate-only* JSON, not original site×species observations.

| Label component | WCVP taxa observed in source | source species×site-years | β for island−mainland light gradient | 95% clustered CI |
|---|---:|---:|---:|---|
| Positive in both | 612 | 8,190 | +0.594958 | [−1.025069, +2.214986] |
| **Old-only positive** | **713** | **8,121** | **−2.496679** | **[−3.678390,−1.314968]** |
| New-only positive | 132 | 1,484 | +0.155299 | [−0.459938,+0.770536] |
| Old published total | 1325 | 16,311 | −1.901721 | [−3.705839,−0.097602] |
| Updated v0.11 total | 744 | 9,674 | +0.750257 | [−1.076060,+2.576575] |

Exactly: −1.901721 = −2.496679 +0.594958; +0.750257 = +0.594958 +0.155299. Thus the new minus old change +2.651978 decomposes into **+2.496679 removal of old-only source labels (94.14%)** and +0.155299 novel source labels (5.86%). Fractions describe model coefficients, not shares of plants or causation.

## Actual old-only plant source versus retained positive plants

| Attribute (same old cache extraction rules) | old-only 713 | retained 612 |
|---|---:|---:|
| Observed WCVP accepted botanical genera | **426** | **399** |
| Source species×site-year units | 8,121 | 8,190 |
| Island species×site-year units | 783 (9.64%) | 570 (6.96%) |
| iNaturalist-titled study records **only** | 402 (56.4%) | 436 (71.2%) |
| At least one old butterfly-positive relationship too | 637 (89.3%) | 534 (87.3%) |
| Old moth relationship labeled `pollinates` | 186 (26.1%) | 137 (22.4%) |
| Median old moth relationship rows / species | 2 | 2 |
| Median independent `study_title` strings / species | 1 | 1 |
| Median documented moth insect names / species | 2 | 2 |

Exploratory Fisher tests using species as rows, three documentation attributes adjusted BH:
- **iNaturalist-only old source:** old-only rate **lower** than retained, Haldane OR≈0.522, p=2.06×10⁻⁸, q=6.18×10⁻⁸. Thus *"the old-only 713 disappeared because all depended on iNaturalist"* is false. The strongest documented class distinction goes in the **opposite direction**.
- Butterfly-positive source claim: OR≈1.224, p=.264 (not resolved).
- `pollinates` database relation: OR≈1.223, p=.124 (not resolved).
- Neither butterfly-positive claim nor GloBI `pollinates` proves effective diurnal compensation, nocturnal pollination, or in-situ pollen transfer.

Old-only species' own source-origin coefficient splits (same original regression and same denominator):
- 402 with **only iNaturalist-titled moth relationships**: β=−1.029236, 95%CI [−1.920913,−0.137559].
- 311 with at least one **non-iNaturalist-titled moth relation**: β=−1.467443, CI [−2.743246,−0.191640].
- Sum −2.496679; both source subsets contribute. These CI/p values are exploratory overlapping-model tests, not independent validations.

## Taxonomic breadth and enrichment testing

**426 genera across the old-only 713 species; 399 genera across retained 612.** Only 143 genera are represented in both sets. Among the 49 genera with ≥5 old-positive taxa across both classes, no post-hoc old-only-versus-retained genus contrast survives BH q<.05. This does not prove equivalence or absence of phylogenetic structure.

Largest old-only genus tallies (old-only/retained):
- *Vicia* 11/3; *Prunus* 10/5; *Trifolium* 10/4;
- *Allium* 9/5; *Rubus* 9/4; *Solidago* 9/4; *Centaurea* 9/4;
- *Vaccinium* 8/0, nominal Fisher p=.0089 but q=.437, not significant after 49-genus screen;
- *Ipomoea* 7/1; *Verbena* 7/1; *Cirsium* 5/10.

These genera span biologically different floral forms and pollination systems, so the data cannot identify a single moth-exclusive syndrome. Genus is a rough taxonomic summary, not a direct flower-color, corolla-shape or flower-time measurement. Correlated species and documentation biases mean Fisher q-values do not establish evolutionary independence.

## Interpretive finding and next proper analysis

The shared **measured commonality** is **being excluded from the updated source-documentation category** combined with a diffuse, negative historical island−mainland association in the old photographic label. The 713 are *not* coherently a single genus, nor demonstrably effective moth-dependent, and the iNaturalist-only hypothesis fails its direct comparator.

To discover an actual biological common trait:
1. Join **WCVP families / genus / branch lengths**, **actual reviewed species floral traits** (flower tube, color, flowering time, day/night visitor effects), and **species-native status on each real island**, with strict data provenance; fd public profile aggregate alone is not enough to label these 713 species.
2. Compare all three positive-documentation groups on a **fixed botanical and geographic study frame** with trait missingness, plant observation popularity and source-study propensity explicitly modeled.
3. Trace the disappeared plant–moth relation **source IDs** through old→v0.11 (same study and taxon IDs) to distinguish removed source edge vs synonym/crosswalk vs changed extraction rule. A later database missing entry does **not** invalidate an earlier observation, let alone prove a biological absence.
4. Independent day/night exclusion and pollen-deposition studies are necessary to claim a moth-dependent ecological mechanism.

**All findings in this note are descriptive/post hoc.** The regression is of source-documented GloBI positivity among photographed plants against same-year light, not local pollinator dynamics or plant decline. Do not replace the preserved README original β or claim a causal ALAN effect from source version turnover.
