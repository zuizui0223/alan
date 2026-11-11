# GloBI-first island/mainland ALAN: real 2019–2025 observational photo-composition test

**2026-11-11 status:** real source observation fit EXECUTED, association is SOURCE-SENSITIVE, and a moth-specific island response is NOT established. New full GloBI v0.11 source scanner is IMPLEMENTED/TESTED on synthetic source but the actual 2.5 GB release has NOT been downloaded or analyzed.

## Which data actually entered this fit?

- The already source-SHA verified 61,967-claim GloBI flower interaction extract (7 moth and 5 butterfly families), NOT the full modern stable GloBI release.
- Original, independently identified WCVP-accepted plant species photographed in same 0.05-degree source cell and observation year during 2019, 2021, 2023, 2024 and 2025. One independent original observation photo represents each observation identity. Observed photos are neither a flora census nor effective moth-pollination observations.
- Same-original-year VIIRS 4-month observation median, WorldPop and bioclimate/elevation sources. DO NOT call the satellite product 12-month annual VNL.
- Outcome: 100 × number of distinct actually photographed accepted botanical species with an independent GloBI-positive moth record / number of all distinct accepted actually photographed species at that photo cell-year. Analogous butterfly output; the difference moth minus butterfly cancels any species that has BOTH positive GloBI groups.

## Geography and model

A source-backed site has at least 3 actual photo opportunities and 2 accepted photographed species (primary source breadth), or 5 photos / 3 species (strict sensitivity). Require an observed island AND mainland site in the exact same observation year and coarse 5-degree geography block. Cluster standard errors by geographic 5-degree block, across repeated years.

Within each year×block, assign equal aggregate weights to observed island and mainland photo cells, avoiding the heavy mainland sampling skew. Fit nightlight, island status, their interaction, log(population density), log(original photo opportunities), log(number of photographed taxa), temperature, log(rainfall) and elevation. Year×block fixed effects absorb time/geography background. Compare using the SAME fixed observed sites and covariates across the four GloBI study-source definitions. Study-source exclusion is a sensitivity, not cherry-picked selection of a positive result.

## Actual numerical result (year×region matched, class-balanced, geographic-cluster t interval)

Main source cohort: 933 island site-years and 9,765 mainland site-years in 78 distinct 5-degree blocks. 195 complete-case year×region groups.

| Source evidence for photographed positive moth plants | Island minus mainland light slope (pp/log1p light) | Nominal 95% cluster CI | nominal p |
| --- | ---: | ---: | ---: |
| All 61,967 GloBI source claims | -1.902 | [-3.706,-0.098] | .0391 |
| Exclude iNaturalist-derived GloBI interaction records | -0.498 | [-1.837,+0.840] | .461 |
| Exclude iNaturalist and trophiCH-derived GloBI records | -0.955 | [-2.245,+0.336] | .145 |

With both positive guild flags on each same distinct plant, the direct **moth-minus-butterfly** outcome has island-minus-mainland nightlight slope -0.778 [-2.701,+1.146], p=.423 (all GloBI). After iNaturalist exclusion the contrast is +0.407 [-0.699,+1.513], p=.466.

Strict photographic cohort: 265 island / 2,658 mainland site-years, 35 independent geography blocks. All-GloBI moth β=-0.930 [-3.161,+1.301], p=.403; iNaturalist-excluded moth β=-0.352 [-2.105,+1.402], p=.686.

**Common nightlight support:** Further require BOTH land classes to have actual overlapping source-light range within year×5-degree block (not extrapolated). Broad cohort reduces to 619 island / 3,468 mainland site-years, 37 blocks. Moth–butterfly contrast: all-source β=-0.532 [-2.895,+1.830], p=.651; iNaturalist-excluded β=+0.959 [-0.071,+1.989], p=.067 (post-hoc, not confirmation).

In 16 balanced models comprising 2 photography gates × 4 positive-evidence source definitions × moth/butterfly proportions, smallest nominal p=.0391 but minimum BH-adjusted q=.626. No adjusted discovery.

**Interpretation:** there is a nominal negative cross-land light-slope difference with the full bounded GloBI sample, but it is not reproduced by source-quality restriction, a stricter photography cohort, or a direct moth-vs-butterfly contrast. These estimates do NOT demonstrate preferential biological loss of actual island moth-pollinated plants, nor a gecko-predation or naturalization mechanism.

## Full GloBI v0.11 2.5 GB source — extraction built, not run

- Zenodo fixed release: https://zenodo.org/records/22691479 (published 2026-09-10).
- Full interpreted interactions.csv.gz SHA256: 670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3.
- scripts/analysis/alan_globi_full_v011_extract.py performs SHA-pinned complete-table chunking, with bidirectional flower-visit/pollination relations; moth is a known lepidopteran family other than butterfly families, and family-unknown remains UNKNOWN. Source study exclusions and insect species-vs-family detection are preserved.
- New workflow .github/workflows/alan-globi-full-v011-manual.yml provides a manually invoked SHA-checked full download/scan. **Not executed in this session.** Manual workflows may require availability on the repository's default branch and suitable Actions permissions; the working Draft PR does not itself imply execution.
- Full GloBI taxa require separate name-resolution/WCVP reconciliation before joining plant photos. New full-version plant species counts and ecological estimates are not currently available. Do not print the old 2,362 positive moth names as a full-v0.11 count.

## Source-specific failures to keep explicit

1. GloBI flower visitation somewhere is not observed visitation at the photographed local island or mainland.
2. A plant without a GloBI moth record is not known to lack moth pollination. The estimated endpoint is documented positive-EVIDENCE fraction, not an unbiased pollinator-guild fraction.
3. Large-scale GloBI meta-source and reference concentration, plant species composition, and differential citizen-science recording can make biological-looking effects appear.
4. Land classes may have nonoverlapping nightlight support. A same-year coarse 5-degree block and a year×block FE do not guarantee the same environmental/cultural baseline.
5. Four-month satellite radiance is not local spectral illumination, neither a year-round exposure nor a treatment assignment.
6. No island-specific gecko predation observations, documented naturalized recruitment dates or reproductive output were joined.
7. Cluster-robust CIs and all p-values are nominal exploratory association diagnostics, not intervention or naturalization inference.

## Reproduce and test

Run the same identity-level source merger documented in the earlier [GloBI-first join protocol](ALAN_GLOBI_FIRST_ISLAND_MOTH_SUPPORT_20261111.md) to produce a **PRIVATE** site-year species-type × original VIIRS table. Never publish that identity/coordinate table.

Then:

    python scripts/analysis/alan_globi_crossland_light_robust_20261111.py --private-source /secure/photo_site_year.csv.gz --out /secure/model_summary.json
    python -m unittest discover -s tests -p 'test_alan_globi_crossland_light_robust_20261111.py' -v
    python -m unittest discover -s tests -p 'test_alan_globi_full_v011_extract.py' -v

The complete model source is stored on ALAN PR #40, with aggregate-only actual outcomes in results/alan_globi_crossland_robust_light_20261111.json. Synthetic unit-test invariants were verified locally. All raw restricted individual source rows stay off the repository. Do not modify fd/island Chapter 1 or overwrite prior ALAN estimates.
