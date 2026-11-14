# ALAN: Does experimental light GEOMETRY reverse reproductive responses? Original Li 2026 audit

## Scientific scope and source lock
**Independent experiment, not the Lu/Tibet nighttime insect visitation dataset.** Li, Adeniji, Meah & Clements, bioRxiv (posted 2026-07-10), DOI [10.64898/2026.07.10.737511](https://doi.org/10.64898/2026.07.10.737511), is **not confirmed peer-reviewed**. Their original `seedset.csv` is available in the [author's publicly readable GitHub repository](https://github.com/Dongboli/experimental-data/tree/main/Artificial%20lights%20at%20night), exact Git blob SHA1 `b616aca5beeba1457549841477d61ca4f796a147`. The author folder identifies this as *data for peer review only*, so no raw observations are redistributed in ALAN; the Actions runner fetches them ephemerally and publishes only aggregates.

Three Bristol/Fenswood hayfields were the original independent spatial blocks, three treatments each: C=no-light; P=isolated low-intensity point lights; L=25 m 26-bulb linear arrays, with positions changed in midseason. Plants of three moth-associated species, *Matthiola longipetala* (ML), *Silene nutans* (SN), *Zaluzianskya capensis* (ZC). **No insect visitors, effective pollen, day/night pollination-exclusion, selfing, predators or plant population light exposure were directly observed.** Treatments also differ in total bulb count and illuminated footprint; array vs point is NOT a geometry-only contrast.

## Original published plant-level data recounted
Original CSV has **79 plant rows**; 78 nonmissing seed-formation outcomes. Original `pollinated` variable is exactly a **binary flag for any `Seeds > 0`**, NOT independently measured pollen transfer. Source true treatment-by-plant records:

| Treatment | Seed-present plants / nonmissing | Seed-presence fraction | Total observed seeds | Total fruits |
|---|---:|---:|---:|---:|
| No-light C | **8/26** | **0.30769** | **672** | **32** |
| Point P | **12/25** (1 unknown) | **0.48000** | **967** | **33** |
| Linear array L | **20/27** | **0.74074** | **3,729** | **149** |

Naive plant-row unadjusted array-minus-dark **risk difference +0.43305**, odds ratio **6.43**, but THESE DO NOT HAVE 79 INDEPENDENT FIELD RANDOMIZATIONS. Point-minus-dark plant-row difference +0.17231, without field/site-robust inference.

## More conservative independence test
Check original three hayfields separately (L vs C, **without pooling individuals across field blocks**):
- Hayfield 1: L **4/8**, C **1/6**, risk difference **+0.33333**.
- Hayfield 2: L **7/9**, C **3/11**, risk difference **+0.50505**.
- Hayfield 3: L **9/10**, C **4/9**, risk difference **+0.45556**.

Thus **3/3 independent hayfields positive**, unweighted field mean difference **+0.43131**, with three leave-one-entire-field-out estimates retaining positive sign. The **minimum two-sided sign fraction for n=3 = 0.25**; it is not an original randomization-inference p-value because the treatment positions were reshuffled mid-season, allocation sequence is unavailable, and unequal flower/plant opportunities persist. Published GLM `p=0.0046` (array versus controls and point model) is the original authors' model output, not a separate inferential claim from this reanalysis.

Check each original plant species separately: ML **10/15 L vs 5/16 C**, SN **5/6 L vs 2/6 C**, ZC **5/6 L vs 1/4 C**, positive source aggregate directions across all three taxa. ZC has no original C observations in hayfield 2; do not infer this missing field-species zero. A secondary entirely present-in-all-fields species subset (ML+SN) retains a positive array-minus-dark observed source difference (**15/21 vs 7/22**, +0.39610). Even these are nested within the same three hayfields; species and plants are NOT independent replicating ALAN experiments.

Original IDs are imperfect: source `Label` can repeat for multiple plant rows; source `id` is missing for most ML plants, and `Label` missing for SN plants. These keys are **not asserted to identify unique independent plants**; do not deduplicate by apparent location label or claim individually randomized 79 experimental units. One P-plant has missing outcome and remains excluded rather than imputed no-seed.

## Relationship to source-original historical ALAN conclusions
The new **positive reproductive direction** is a deliberately scoped counterexample to a universal negative ALAN impact:
- Knop et al. (2017), [Nature DOI 10.1038/nature23288](https://doi.org/10.1038/nature23288): original Swiss field manipulation reported **62% reduction in nocturnal visits**, **13% reduction in one focal plant's fruit set** with diurnal visits present. This is the ORIGINAL authors' result; a new matched seed-data reanalysis of this study is not supplied here.
- Lu et al. (2026), [Ecology and Evolution DOI 10.1002/ece3.74368](https://doi.org/10.1002/ece3.74368): original Tibetan five-field study, **57.4% fewer nocturnal Lepidoptera visits**, more Diptera but **no directly measured fruit/seed outcome**. Source-original independently reproduced in ALAN PR #65 and cross-sheet discrepancy flagged in PR #72.
- Boom et al. (2020), [Scientific Reports DOI 10.1038/s41598-020-75471-1](https://doi.org/10.1038/s41598-020-75471-1): original Netherlands *Silene–Hadena* nursery-pollination spectral treatment; original study's overall estimated seeds/flower test **p=0.63** despite some component changes, ALAN PR #67.
- Li et al. (2026), **this actual original 79-plant source**, 3 spatial fields: array reproductive response positive; isolated point effect less resolved. Night-visitor behavior was NOT recorded.

Because these four studies use different **visitation vs fruit set vs seed count** endpoints, plant syndromes, lighting spectrum/geometry, and units, there is no valid four-study common-effect meta-analysis on a single nightlight coefficient. The testable macro hypothesis shifts to a **light context × realized functional pollinator dependency × resource opportunity** model, with direct reproductive outcomes, not generic satellite brightness as a sole ecological mechanism. The hypothesis that arrays act as movement corridors is attractive but **not demonstrated** by the author's seed table.

## Real-source reproducibility
The executable test downloads only the original author CSV at runtime and SHA1-validates its exact Git blob, verifies 79 source rows, 78 measured outcomes, three hayfields and three species. It computes arm, field, species, round and missingness-aware contrasts with synthetic tests. See `scripts/analysis/alan_Li2026_real_geometry_seed_hayfield_source_20261114.py` and `.github/workflows/alan-Li2026-genuine-geometry-seed-threefield-20261114.yml`. Its uploaded result is only a machine-readable aggregate and original-source SHA receipts; no raw source coordinates/plant labels are published.
