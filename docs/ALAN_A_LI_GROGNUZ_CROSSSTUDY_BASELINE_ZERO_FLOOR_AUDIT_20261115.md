# ALAN A cross-study reproductive hurdle hypothesis: two original light experiments, baseline-floor falsification (2026-11-15)

## The new, source-proven result

**A dark control with measured zero seeds cannot have a negative light-minus-dark difference when all realized seed counts are nonnegative.** Hence "light improved fertility in dark-zero strata" is *arithmetically selected*, not in itself a biological rescue effect. The larger ecological claim that ALAN systematically benefits poorly reproduced flowers **does not survive as an established result**.

We independently reran two actual source datasets with appropriate scope: Li et al. 2026 moth-syndrome **preprint** with open nighttime flower exposure (original 79 plants), and Grognuz & Knop 2026 **preprint** with flowers **bagged during three prior illuminated nights and then offered 90 minutes daytime-only visitation** (source-authenticated supplementary Table S3 block means). Both are independent **papers**, but their species, light sites, flower access, denominators, measured outcomes and assignment blocks differ. No new pooled meta coefficient, p value, island interaction or nocturnal pollination mediation is estimated.

## Li 2026: three-arm control, not just light vs dark

Original author Git blob `b616aca5beeba1457549841477d61ca4f796a147`, same field × species × flowering round. Exactly **11** strata with **all three measured treatments**: L=25m linear array (26 lamps), P=isolated point lights (4 lamps), C=dark. 5 of the 11 have zero measured dark seeds.

| Original within-stratum comparison | Pooled seeds/starting bud RR | Signed contrast count | True independent field mean signs |
|---|---:|---|---|
| **L/C** | **3.426** | 9 positive, 2 negative | 3 positive, 0 negative |
| **P/C** | **1.052** | 4 positive, 4 negative, 3 ties | 1 positive, 2 negative |
| **L/P** | **3.256** | **8 positive, 3 negative** | **2 positive, 1 negative** |

The 3-arm positive L/P direction **is not entirely forced by dark zero** because C is absent from L−P. It *still cannot establish a geometry-only light effect*: P and L had different lamp numbers, footprints and light fields, and the original illumination unit was only **three independent meadows**. Restricting to **at least two plants in EACH L/P/C source stratum leaves 6 strata**, with seed/bud L/P **3 positive/3 negative**. The earlier 12 C/L stratum analysis had 6 dark-zero strata; the current 11 matched **three-arm** stratum subset has 5. Do not confuse the two selection sets.

**Li is a genuine positive reproductive experiment under a linear array but not evidence that all ALAN benefits plants, rescue occurs at dark-zero sites, or moth pollen movement is the mediator.** The paper itself calls for direct moth tracking. There is no correct site-randomized CI for the L/P contrast recovered here.

## Grognuz & Knop 2026: a different reproductive pathway

The public supplement Table S3 used here contains four ALAN × slug treatments per daytime choice block, summarized as source-observed per-four-pot group means and observed denominators, not authentic individual pot-level raw rows. **Twenty Trifolium pratense blocks** and **eight Centaurea jacea blocks** reuse light installations and some plants. True lamp assignment clusters are not 28.

| Plant | Choice blocks | ALAN−dark seed count/target mean difference (slug arms equal weight) | ALAN positive blocks | ALAN negative blocks | Dark-zero blocks | Excluding dark-zero, ALAN−dark mean (posthoc, outcome-selected) |
|---|---:|---:|---:|---:|---:|---:|
| *Trifolium pratense* | 20 | **-1.406** | 8 | 12 | 1 | **-1.564** |
| *Centaurea jacea* | 8 | **-0.422** | 3 | 5 | 1 | **-0.518** |

Both dark-zero blocks were positive by arithmetic nonnegativity; the overall plant-specific means were **negative** and became more negative when the dark-zero block was omitted. This **is not a valid outcome-adjusted causal estimator**: excluding strata after observing control seeds causes selection bias, and the published Table S3 blocks are dependent through lamps and pot reuse. The original authors' statistical model, published light effect and missing outcomes must remain distinct from these signed means.

## Why the naive "dark baseline × light effect" finding is dangerous

Define C = independently observed control outcome and L = illuminated outcome. For the **observed difference** D = L−C,

`Cov(C,D)=Cov(C,L)−Var(C)`.

Even under a **true null** with independent equal-distribution L and C, this is `−Var(C)`, a negative association. More strongly, if `C=0` and all seed counts are nonnegative, `L−C≥0` is a logical necessity. These patterns need NOT be caused by pollination rescue or an underlying trait moderator. Negative covariance from the shared observed control is similarly relevant to the primary regression combining nighttime V and lit/dark y.

Therefore future A meta regression **must use pre-treatment independent control fecundity or independent access experiments**, or explicitly model shared-control sampling covariance; post-treatment observed dark counts cannot be used as a naive explanatory covariate for the same light contrast.

## Frozen meta-analysis disposition

* **P3 total reproductive ALAN pooled effect:** not estimated. Li seeds / **initial buds** vs Grognuz seeds / **target flower or inflorescence** under daytime-only pollen access cannot be treated as one estimand.
* **Main A outcome × functional nocturnal dependence V meta regression:** 0 confirmation-ready studies in currently audited nine directly manipulated plant reproduction publications. Lu open seed response and covariance unavailable; Macgregor original day/open cell data unavailable.
* **Island×ALAN reproductive interaction:** no matched experimental replication here. Trait-based global GloBI/flower photographs remain observational.
* **Recommended next empirical recovery:** true lit/dark × (day/night/open/closed) plot-pair original outcomes for Lu and Macgregor, and first-flower marked mature fruits/total seeds plus cluster identifiers for Knop, Boom, Giavi, Li, Grognuz. Keep published positive, negative, zero and unknown directions; no p-value selection.

This is a **posthoc, source-audited evidence map and falsification exercise**, not a registered meta analysis. It establishes a specific non-identifiability to avoid in further global ecological theory; no original author raw source was redistributed.

## Reproducibility

[Exact source-check Actions #43722481561](https://github.com/zuizui0223/alan/actions/runs/43722481561) succeeded with **4/4 tests**. It calls [Li SHA1-verified source triarm code](../scripts/analysis/alan_Li2026_original_triple_arm_20261115.py) and [Grognuz original-derived supplementary block guard](../scripts/analysis/alan_Grognuz2026_block_floor_signs_20261115.py), then the [combined nine-paper eligibility guard](../scripts/analysis/alan_A_crossstudy_baseline_floor_gate_20261115.py). [Machine-readable signed summary](../results/alan_A_crossstudy_hurdle_source_signed_summary_20261115.json); full result artifact retained in the workflow run, no author plant/site raw data committed. Main remains unchanged; same PR #78 Draft.
