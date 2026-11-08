# ALAN — island area × nightlight: photographed taxon assembly vs taxon-relative ROI colour

**Source-locked analysis protocol, 2026-11-08.** New empirical estimates will be appended only after the actual ALAN workflow succeeds.

## Research question

The previous real-data study found a positive island-area × nightlight interaction for observed white-colour mass among photographed floral regions, after controlling for mainland shoreline isolation. Why can the *photographed* floral-colour assemblage differ with area and nightlight?

1. **Photographed taxon assembly:** distinct WCVP accepted plant taxa appear in the verified original-photo sample of different island sites. Use the taxa's photo-archive mean nine-palette white fractions to describe the component corresponding to photographed accepted-taxon identities.
2. **Taxon-relative photographic variation:** among the actual photos of a given accepted taxon, the measured flower-region white fraction differs from its photographically observed reference mean. This includes potential within-species heterogeneity **and** differences of lighting, flower-part ROI, phenology and volunteer photo choice.

A species' **photographic mean colour is not an independently verified invariant botanical trait**; an archived mean derived from the same photos is an algebraic reference, not proof of field community turnover or adaptation. Nondetection of a taxon in images is not field absence.

## Fixed source and original comparison

- The exact SHA-verified original 2019, 2021, 2023, 2024 and 2025 photographic ROI9 and WCVP accepted taxa, coupled to genuine year-specific VIIRS and WorldPop.
- The exact same GSHHG island selection as the earlier area/isolation study: **44 islands / 650 qualifying photo cell-years** in the strict 5-photo/3-accepted-species cohort and **95 islands / 1,884** in the broader 3-photo/2-accepted-species cohort.
- GSHHG area and distance to mainland coast vertices, with **island and observation year fixed effects**, population, BIO5, rainfall, elevation and original photo effort covariates, and island-cluster nominal uncertainty.
- Two prespecified weights: equal photo cell-year, and island-equal weights summing to one per GSHHG island.
- **Exact coefficient identity**: with the same rows, design matrix and weights, both interaction estimates for measured white ROI must equal the sum of its taxon-reference component and taxon-relative photographic deviation component. Baseline signed area and mainland-distance interaction coefficients must match previously published ALAN numbers within numerical tolerance.

## Independent-year diagnostic (different target population)

- **Training photo years:** 2019, 2021, 2023.
- **Evaluation photo years:** 2024, 2025.
- Only genuinely observed accepted taxa having prior-year photo colour references enter the 2024/25 comparison. Require at least **two** earlier-referenced accepted taxa per measured test photo cell and at least **two distinct original 0.05° sites** on an included GSHHG island.
- Evaluate the 2024/25 observed colour of only those matched taxa, prior-year species-colour photographic reference, and their exact difference; use the same island×area and island×distance modeling form if **at least 20 independently supported islands and 100 eligible test photo cell-years** remain. Otherwise record actual support rather than fitting an implausible model.
- Unlike the full-source algebraic decomposition, future photos do **not** define their own reference. But the subset of taxa with records in earlier years is selected and need not represent the whole flora.

## Interpretive ceiling

If the archive-based taxon-reference component covaries with island size and nightlight, this supports a **photographed taxon composition-related association**, not a measured ecological species replacement mechanism. If taxon-relative photographic colour does, this is not automatically genetic flower-colour evolution or nocturnal pollinator selection. VIIRS upward radiance does not measure ground-level flower illumination.

This study does not use speculative missing-colour corrections or simulate new species. Photographs, original colours, geographic names and historic results remain unchanged.

## Reproduction entrypoints

- Code: `scripts/analysis/alan_island_area_taxon_photo_decomposition_20261108.py`
- Tests: `tests/test_alan_island_area_taxon_photo_decomposition_20261108.py`
- Workflow: `.github/workflows/alan-island-area-taxon-composition.yml`
- Source base: [ALAN combined island area versus mainland isolation](ALAN_MAINLAND_ISOLATION_VS_ISLAND_AREA_EMPIRICAL_20261108.md).

The empirical coefficients and signed action artifact are not asserted until the workflow reports success.
