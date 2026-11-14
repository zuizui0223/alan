# ALAN — genuine independent Swiss and Tibetan guild responses after source-measured flowering cover

**New cross-source method, 2026-11-14.** [Original field data QA and final same-source rerun SUCCESS #43574048558](https://github.com/zuizui0223/alan/actions/runs/43574048558), [aggregate-only result #14105441985](https://github.com/zuizui0223/alan/actions/runs/43574048558/artifacts/14105441985). This uses two **independent original artificial-light experiments**, not two independent global biomes sampled in a common design. Swiss original Giavi, Fontaine & Knop (2021), *Nature Communications*, [DOI 10.1038/s41467-021-22011-8](https://doi.org/10.1038/s41467-021-22011-8), authenticated Zenodo [4540407](https://zenodo.org/records/4540407); Tibet Lu et al. (2026), *Ecology and Evolution*, [DOI 10.1002/ece3.74368](https://doi.org/10.1002/ece3.74368), authenticated nested Europe PMC workbook (actual SHA256 ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f).

The **original real source reanalysis** independently reconstructs: Swiss **6** paired lit/dark field sites, **37** paired survey occasions, 21 plants, 2,384 actual documented interactions in 1,884 source rows; Tibet **5** independent original blocks, four rounds, 400 true-day/night group-count records (night and day each have 100 lit and 100 unlit). Plant taxa are not assumed the same between regions and study-site nested records do not become 1,884 independent sites or 400 randomized plot replicates.

## Main genuine original-source result

**Swiss possesses measured plant floral COVER proxy (plab), Tibet does not.** We deduplicate Swiss plab to once per original plant × paired sampling event × light treatment × site, instead of treating its three guild-specific copies as independent flower measurements. There are **628 original plant-cover source records** and **zero within-same-plant/event disagreements** after cross-guild reconciliation. The proxy is a discrete estimate of floral cover, *not* direct nectar availability, actual flower opening time or same-day pollinator activity overlap.

For a given guild g at each independent field site or block, calculate:

- **Unadjusted source-visitation rate ratio**: (Lit visits + 0.5)/(Unlit visits + 0.5) × (Unlit field survey opportunities/Lit opportunities). 0.5 is a continuity correction for zero counts, so a count-exact 2:1 can be slightly below 2.
- **Swiss floral cover-adjusted ratio**: same count ratio × (measured Unlit plab total / measured Lit plab total). The denominators are source **observed floral-cover proxy**, NOT a true census of every open pollen-receptive flower.
- Average **log ratios equally over the six Swiss independent sites or five Tibetan independent blocks**, exponentiate to a geometric ratio. Uncertainty: 5,000 **whole-independent-site/block resamples**, descriptive quantiles only. Do not use site×plant×insect rows as independent effects.

### Two-source original daytime guild ratio comparison

| Real original daytime guild | Swiss, 6 sites, visits per paired survey | Swiss, 6 sites, after observed flower-cover adjustment | Tibet, 5 blocks, visits per original survey |
|---|---:|---:|---:|
| **Diptera** | **0.971** | **1.018** | **1.471** |
| **Hymenoptera** | 0.900 | 0.943 | 0.792 |
| **Coleoptera** | 1.191 | 1.248 | 0.377 |

Numbers above are independent-site equal-weight, original source group-count ratios, not fitted cross-study pollen or actual fruit/seed causal effects. The original authors' pooled raw group counts can give different signs because a single high-count site dominates. Swiss raw **Diptera** 569 dark →550 lit **−3.34%** (pooled), yet 5/6 sites have a positive *raw* within-site change; after deduplicated floral cover standardization the ratio is 1.018 and only **3/6 sites are positive** (3 negative). Its six-field bootstrap descriptive 95% CI for the adjusted ratio is **[0.419, 2.317]**, spanning strong loss and gain. Swiss Hymenoptera adjusted ratio CI [0.560, 1.463] and Coleoptera [0.423,3.049] similarly unresolved. The geometric mean independently observed Swiss *lit/dark floral cover* ratio is approximately **0.954** (the shared adjustment factor in the three guild results). This is source availability/composition, **not** a measured causal effect of ALAN on Swiss flowering onset.

Tibet daytime Diptera independent-block equal-weight ratio is **1.471**, 5/5 blocks show a higher lit visit rate, and descriptive block bootstrap 95% interval **[1.288,1.756]**. Hymenoptera 0.792 [0.689,0.907], Coleoptera 0.377 [0.223,0.638], each 5/5 source blocks lower. These are robust internal descriptive block trends from an already published experimental light treatment, **not** five independent publications and not a global meta p value. A two-sided complete five-block sign test has p ≥ .0625 at best; the source's original mixed-model effects should remain the formal original inferential source.

### Real guild balance and composition (same three insect orders)

- Tibet daytime **fly/bee relative visitation-channel balance** (ratio of lit/dark ratios) rises **1.858×**, **5/5** blocks. Swiss equivalent **1.079×**, **4/6** sites, extremely heterogeneous site-specific ratio range 0.103–3.200. These group differences *cancel* a shared floral-cover denominator, but remain visit-count ratios rather than different insects' conspecific pollen deposition efficiencies.
- Same exactly three orders Diptera/Hymenoptera/Coleoptera: original independent site-mean Jensen–Shannon **guild visitation composition** divergence under lit versus dark is **0.0613 bits** in Swiss day, versus **0.0327 bits** in Tibet day. Their identical taxa make this *numerically* comparable as descriptive compositional statistics, but two ecosystems have non-matched species pools, light spectra and survey geometry; no ecological effect ranking or inferential cross-region means are claimed.
- Tibet night may have substantially different 3-order composition change, but **the three shared orders EXCLUDE dominant nocturnal Lepidoptera**. This omission means a day-versus-night contrast of the *three-order* JS indices does NOT capture the full pollination network. A separate fully documented Tibet **five-order** day/night original within-block diagnostic is calculated in this branch, without treating actual light–flower phenology as measured.

## Original Tibetan complete five-order *night-versus-day* source result

For Tibetan original five field experimental blocks only, compute treatment-versus-control JS divergence on **ALL FIVE** source order columns, including the nocturnally dominant Lepidoptera. Compare the two diel periods within each identical experimental block, not across differing ecosystems or Swiss plant pairs.

| Experimental diel period | Five-order independent field-block mean JS divergence (bits) |
|---|---:|
| Diurnal | 0.04483924 |
| Nocturnal | **0.07298861** |
| Original mean of within-block nocturnal minus diurnal | **+0.02814937** |

Nocturnal light/dark guild-composition divergence was larger than diurnal in **5 of the 5** actual blocks (paired difference range +0.01886 to +0.04236 bits). This constitutes an *observed five-block compositional pattern within ONE original study*. Its minimum two-sided sign fraction is 2/32=.0625, and it must not become a new significant global meta-analytic discovery.

**Crucial biological interpretation:** This difference is a taxonomy-specific change in recorded interaction composition across two 4-round daily survey windows. It does not measure functional reproductive pollen delivery, existing island moth pollinator shifts, predators, floral opening/closing, seasonal flowering–visitor overlap, or plant seed formation. Main total vs five-order counts differ slightly in 23/400 original Tibetan rows; five-order JS is computed on explicitly recorded guild counts without silently reallocating the discrepancy.

## Why this changes the global syndrome question

The two source experiments do not support an inherent rule that "flies can rescue moth pollination under ALAN" simply because fly visit counts rise in one alpine meadow. Within the other Swiss experiment, even the **direction of fly visit-rate response per measured flowering resource is not stable across six field pairs**, despite the known positive five-of-six *raw* site count direction. Plant abundance and independent geographical opportunity are critical even before effective pollen movement, predators or seed counts enter the model.

New working target for a truly MACRO analysis is therefore:
1. Measured **per-open-flower effective conspecific pollen delivery** for nocturnal Lepidoptera, diurnal bees and other vector guilds under paired ALAN exposure;
2. Source-measured plant breeding mode and autonomous selfing, and independently measured seed/fruit reproduction, including **nursery pollen-vector/larval seed-predator tradeoffs**;
3. Season/day flowering phenophases and pollinator activity with shared site-year and equal time-bin detection opportunities; distinguish *floral availability* from iNaturalist observer DOY;
4. Existing island pollinator-shift channels (moth, bird, bat, gecko, honeybee, fly) actually observed locally before illumination, compared with phylogenetically/regionally matched mainland floral communities;
5. Direct adult-moth predation risk, florivore removal and the predator's actual **positive pollen deposition** as separate mediator paths, not GloBI dual-role species-name labels.

Until these conditions are met, the original 562 islands / 14,330 flora photographs and the independently original 301 diel reproductive field outcomes cannot be treated as a single site×year×taxon **functional interaction** regression. This branch DOES NOT produce a global ALAN×syndrome×predator×phenology×island coefficient, because the required common-support evidence is genuinely absent.

## Reproducibility

Source script: scripts/analysis/alan_two_region_original_guild_flowercover_20261114.py; tests: tests/test_alan_two_region_original_guild_flowercover_20261114.py; original GitHub Actions workflow: .github/workflows/alan-two-original-regions-guild-flowercover-20261114.yml. Source ZIP/DAT files are fetched at run time and authenticated against publishers' fixed original SHA/MD5. No individual animal×plant×source-site record is copied into the public artifact. All previous PRs and main branch remain unchanged.
