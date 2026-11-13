# ALAN 39-island day/night functional phenotype: source stability and direct selfing

**2026-11-13: real original-source rerun successful.** Draft PR #47. Earlier PR #43-46 and main unchanged.

- [Initial exact-source baseline + independent reproductive traits: run #43365293498](https://github.com/zuizui0223/alan/actions/runs/43365293498)
- [Additional genus, continuous, effort weighting and hours sensitivity: run #43365792660](https://github.com/zuizui0223/alan/actions/runs/43365792660)
- [Final aggregate-only JSON artifact #14010366538](https://github.com/zuizui0223/alan/actions/runs/43365792660/artifacts/14010366538)

## Original observational source population and inferential ceiling

- Original 93,209,351 photographic observation IMAGE IDs, not 93m independently analyzed pollen-carrying flowers; 16 signed WCVP photo-flora SQLite groups.
- Historical Li–Cao 1992–2024 PCNL original series, annual mean log1p exposure only up to the photo observation year minus 1, verified source SHA.
- Exact source: 562 genuine GSHHG islands, 77,647 photographed 0.05 degree site×years, 14,330 WCVP accepted botanical taxa, 2,012,925 species×site×years. Original photographed years: 2019, 2021, 2023, 2024, 2025.
- External Kendall & Nicholson 2025 original 1,652 day/night treatment records, 262 original matched study×species×outcome pairs and 228 comparable fruit/seed-setting paired metrics from 131 independent studies. Species-level daytime/nighttime numerical advantage is a categorical proxy for externally observed reproduction, not local moth pollination.
- Study matched 52 botanical taxa occur in the original photographed-island cohort. The source is not a botanical census of wild plants, and absent study labels do NOT prove biological absence.

## Main source model: exact common-pollination-documentation support

The original post-hoc sample admits photographed site×years with both an experimentally nighttime-higher AND daytime-higher floral species, and at least two geographically distinct light values within the same island×year.

This yields **2,373 photographed cell×years from 39 independent islands**, β = **−1.1747574942** pp per prior mean log1p PCNL, 95% island-cluster CI **[−2.47017,+0.12066]**, nominal p **.07422**. Island×year fixed effects, log photo effort, equal island total weights. Not statistically resolved.

### Exhaustive island, species and genus omission (never choose only favorable omissions)

| Sensitivity | All trials | Estimable | Negative coefficient | β range pp | Nominal 95% CIs fully negative |
|---|---:|---:|---:|---|---:|
| Leave one island out | 39 | 39 | 39 | −1.4285 to −0.6947 | 1 |
| Leave one experimental plant species out | 52 | 52 | 52 | −1.8650 to −0.4874 | 3 |
| Leave one botanical genus out | 42 | 41 | 41 | −1.2732 to −0.5479 | 1 |

One whole-genus deletion leaves <30 independent islands, so no false nominal p-value is reported. Many source plant deletions change which originally photographed cell-years remain eligible. These overlapping source deletions do NOT amount to independent confirmatory experiments or proof that the true ecological effect is negative.

### Stronger alternate measurement checks, crucial to scientific interpretation

**On exactly the same original 39 islands / 2,373 site×years:**
- Use original **continuous author-experiment nighttime reproductive share** instead of the dichotomized night/day label: **β +0.06021 [−0.48877,+0.60919], p=.8255**. The binary negative association is not reproduced by the continuous measure.
- Keep the original binary outcome but give **each photo site-year equal weight** rather than each island equal total: **β −0.22098 [−0.53358,+0.09162], p=.1606**. Substantial attenuation.
- Define species' night/day dominance after correction for actual recorded **treatment hours**: 39 photographed experimental taxa have the necessary source data. A revised co-documented cohort has 436 raw paired-class cell-years in 51 islands, but only **343 truly same-island-year light-qualified cell-years on 12 independent islands**, **NOT ESTIMABLE** at the locked >=30 independent island rule.

**Adjudication:** the sign of a categorical photographic-composition coefficient is stable to one-at-a-time island, species, and genus removal. It is NOT robustly replicated when replacing that label by the actual continuous reproductive contribution, and the hours-qualified cohort is too sparse. Global ecological effect remains unresolved.

## Independent direct selfing and mating systems

The independent 2026-09-08 Island author-reviewed, species-DIRECT high/medium original reproductive-trait ledger was SHA-verified (134,268 source trait rows/62,423 species across all axes; checksum 7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577). Each trait is admitted only as a species-direct single valid value. Distinguish SC from autonomous selfing.

| Scope | Autonomous selfing present/absent | SC/SI | Predominant mating system |
|---|---|---|---|
| All 134 independent externally studied reproductive species | 13/21; Fisher OR 2.13, p=.481 | 29/22; OR .813, p=.782 | Only 1 selfing vs 16 outcrossing; not estimable |
| **52 species photographed in real 562-island cohort** | **10/13; OR 2.40, p=.414** | **22/13; OR .972, p=1.000** | 0 predominantly selfing vs 14 outcrossing; not estimable |

The smaller 52-study-species frame has **23** eligible autonomous-selfing positive/negative taxa, **35** valid SC/SI taxa and **17** mating-system names (3 mixed), far below support for a geographical ALAN × selfing effect. Recorded SC does not demonstrate reproductive assurance. Lack of association is NOT demonstrated equivalence or biological independence.

## Scientific decision and next data gate

Preserve descriptive positive result (stable categorical point sign across one-at-a-time perturbations), but do not promote the exploratory p=.074 to a claim of ALAN-induced moth pollination collapse, population disappearance or compensatory selfing. Source selection, plant genus, native/cultivated status, actual pollinator and reproductive field measurements, and variation in historical light product remain.

Future work should prioritize independently observed **effective day/night pollen deposition, fruit/seed response, autogamy and site-level dark refuge** for the SAME native populations on replicated islands; no assumption that more flower photos by themselves can prove nocturnal effective pollination.

Source code, unit tests and signed source-replay workflow are in Draft PR #47. Result artifact exports aggregates ONLY, not source restricted plant trait–island/coordinate identities.
