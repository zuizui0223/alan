# ALAN — original direct population-year visitor × fruit-set bridge (14 November 2026)

**REAL SOURCE SUCCESS**: [GitHub Actions #43560918822](https://github.com/zuizui0223/alan/actions/runs/43560918822), [original aggregate artifact #14099294494](https://github.com/zuizui0223/alan/actions/runs/43560918822/artifacts/14099294494), [Draft PR #62](https://github.com/zuizui0223/alan/pull/62). Both unit tests and genuine original-source jobs succeeded.

## Real public biological source

This audit goes beyond checking whether a paper mentions a visitor. We connected independent raw ORIGINAL plant-level exclusion outcomes and ORIGINAL visitor counts from Wenzell et al. (2024, *Ecology and Evolution*), [DOI 10.1002/ece3.11443](https://doi.org/10.1002/ece3.11443), [public author data](https://github.com/KWenzell/Castilleja_sessiliflora_exclusion) to the 301 authentic study/plant/actual field site-year fruit-seed pairs of Kendall & Nicholson (2025, *Ecology Letters*, DOI 10.1111/ele.70036).

The global meta-analysis original author treatment file remains SHA-verified as 1,652 treatments / 140 original study IDs, Git blob 39a82d966afcdde670603add1c881ace8791ffe9. The pre-light original 301 paired field source is SHA256 21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af. Independent Wenzell published author original files and actual original Git blob checks are frozen:
- plant-level fruit/flower treatments: 85491786cd9381120bea598b9a9d73ae08acd5ad
- true population/year/visitor-majority labels: 6a620db6f3ffb4425ff48b114056ecf9907dbb5e
- observed flower-visitor counts in two separately collected narrow/wide protocols: 050f2db829a7ce06c7007d85fbb8c4c1010122d1

No post-experiment nightlight year or pre-published manuscript date is substituted. This biological source has NOT experimentally manipulated nightlight.

## Verified original source cross-walk

**Wenzell 2019 direct field experiment:** 7 distinct year×population combinations from 3 publicly studied Castilleja plant taxa, **204 original actual single-plant treatment rows** across day-only, night-only and open treatment arms. Original visitor event records identify majority activity for each population-year: 3 nocturnal hawkmoth-dominant and 4 diurnal bee/hummingbird-dominant. A total of **750** original flower-visitor events appear across the separately collected source narrow + wide visitor observation records, with source population-level totals verified.

**Global Kendall/Nicholson 301-outcome source contains six of these seven true experiment population-years**, all labelled 2019 and exactly corresponding to the original botanical plant name. After independently joining author-reported coordinate, accepted source plant and year:
- 5/6 mapped populations are within **1 km**; conservative strict source tier.
- 1/6 is about **1.65 km** away; **provisional 1–2 km tier**, not claimed an exact coordinate match.
- 6/6 original global study outcome signs for DAY vs NIGHT reproduced using the independently archived direct 2019 plant-level outcome data, despite **not-identical aggregation units** (global author original reported treatment means vs new sum of fruits divided by sum of flowers).
- Maximum source-global versus direct night/(day+night) fraction difference is **0.04558**, so do not claim their numerical means are identical.
- The seventh direct original 2019 population is genuine published information absent from the frozen 301, **not** an extra newly independent journal study or a silently appended data row.

## Original descriptive ecology with authentic same-population/year animal guild

Night reproduction share is calculated for each population-year as actual NIGHT-only fruit/flower success divided by the sum of NIGHT-only + DAY-only measured fruit/flower success. Each published population contributes equal weight, not each plant or visit count.

| Source original comparable set | Nocturnal dominant source populations | Nocturnal mean NIGHT fraction | Diurnal dominant source populations | Diurnal mean NIGHT fraction | Descriptive fraction contrast (night majority−day majority) |
|---|---:|---:|---:|---:|---:|
| All seven independently published 2019 populations | 3 | 0.69845 | 4 | 0.24941 | **+0.44904** |
| Six actually in original global 301, <=2 km (one provisional) | 3 | 0.69845 | 3 | 0.25999 | **+0.43846** |
| Five strict source/global geographic joins, <=1 km | 3 | 0.69845 | 2 | 0.20888 | **+0.48957** |

For the 7 true original field population-years, **6/7** show agreement between the independently observed predominant visitor activity period and the day/night treatment with greater measured fruit set. **One hawkmoth-dominated short-tubed population** nonetheless has higher measured daytime fruit set, a real exception already described in the original author article. For the 6 global source matched source populations, directional agreement is **5/6**, whereas *global source vs direct data treatment contrast sign* matches **6/6**.

The globally frozen original source day/night treatment fraction, grouping the six matched populations by their independently documented source visitor majority, is 0.69480 for hawkmoth-dominant and 0.28766 for diurnal-dominant populations (descriptive across the ONE original independent paper). Thus the direct original data and the global meta-analysis author treatment summaries are directionally concordant without being numerically identical.

## Exactly what is and is not verified

Verified:
- a **population/year** visitor guild predominance measurement exists and matches direct fruit/seed fitness treatment outcomes for one original source study,
- original authentic global/primary author study year/species/geography is 2019, not merely GloBI species-positive at some unrelated place,
- original source experiment's specific NIGHT vs DAY fruit-set sign agrees for six independent source population joins,
- strictly original 1 km geography supports five, while a 1–2 km candidate is separately flagged.

NOT verified:
- the recorded visit occurred on the exact bagged individual plant, or exactly on the 48-hour intervention dates;
- a positive visit necessarily deposited pollen on a stigma or caused fruit set;
- the six or seven field populations are separate original independent published studies — **they are one publication**, and cannot support a cross-study p value;
- other global 301 source papers contain comparable original matched visitor data until separately retrieved and verified;
- anything about nightlight reducing moth, flower evolution, replacement service under ALAN, species loss or island-specific pollination.

**Scientific novelty and scope:** Wenzell et al. already published this hawkmoth-versus-secondary-diurnal-pollinator reproduction function. The incremental result here is the independently **source-authenticated join to the global 301 genuine field reproduction framework**, including the exact number of matched and nonmatched source populations and the geographic source error tier. Do not present the underlying ecological contrast as a new discovery or infer an ALAN effect.

## Future source expansion

This is the first confirmed source-study **population-year experimental reproductive outcome × measured local visitor guild** bridge in the existing 301 empirical sample. The other 91 flagged observer-publications remain unverified at this exact study/time grain. Next attempts should prioritize primary data repos or supplemental visitor tables and require the SAME original source study, accepted plant name, population site, original experiment year, then distinguish actual visit, pollen transfer, and measured fecundity. Require >=25 independent original studies and >=20 plants for any confirmatory guild-comparison inference.

**Authoring:** scripts/analysis/alan_real_Wenzell2019_populationyear_pollinator_reproduction_20261114.py, tests/test_alan_real_Wenzell2019_populationyear_pollinator_reproduction_20261114.py and .github/workflows/alan-authentic-Wenzell2019-guild-reproduction-20261114.yml. Aggregate-only numerical result in original Actions artifact #14099294494; no plant×source site×coordinate rows exported.
