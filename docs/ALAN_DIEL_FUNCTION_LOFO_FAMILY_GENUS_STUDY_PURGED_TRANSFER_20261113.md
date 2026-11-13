# ALAN: Can night-flowering reproductive function transfer across botanical families?

**Original source run SUCCESS on 2026-11-13.** [FINAL original-source GitHub Actions 43403060451](https://github.com/zuizui0223/alan/actions/runs/43403060451), [aggregate JSON artifact 14025518642](https://github.com/zuizui0223/alan/actions/runs/43403060451/artifacts/14025518642), Draft PR 55.

## Original data and inferential ceiling

Original real paired night-only versus day-only pollinator exclusion experiments (Kendall & Nicholson, Ecology Letters 2025, DOI 10.1111/ele.70036) provided 301 true source plant×study×site×year×fruit/seed outcome pairs, with 118 independent published studies and 123 botanical taxa. Independently pinned 1992–2021 original Li–Cao annual satellite PCNL were sourced from Zenodo with all annual MD5 verified and strictly previous to field year. Author 142-species independently compiled floral trait file blob SHA 18381e550a7e405fe4f198909f0bb381bbecf1e0. Night-only reproduction is NOT a direct identification of moth pollinators; satellite light is NOT manipulated lux at flowers.

All four source feature models use the **identical 196 original reproductive site/year pairs with complete documented night/day bloom time, scent, nectar, recorded color, flower shape and positive flower length**, representing **73 studies, 74 plant species, 28 families and 48 genera**. Missing bloom time never assigned day. Original previously verified study-year/geo/elevation/metric-adjusted night-bloom versus day-bloom actual reproductive share: **+0.2147499065 fraction** (+21.475 percentage points), source study-cluster CI [+0.1075,+0.3220], independently reproduced before new testing.

Four originally preselected models: source geography/year/PCNL/elevation/fruit-or-seed metric only; context plus night blooming; context plus scent/nectar/white/tubular shape/log flower length; context plus all floral traits. Train-only standardized fixed Ridge alpha 5; source publication-equal weights recomputed WITHIN each training fold after full family/genus and held-out-species removal, nonnegative output clipped to genuine night-share [0,1], no retuning on held-out families. No original study, species or family identity used as predictor.

## 28 leave-one-whole-botanical-family-out predictions

Each source family is fully withheld exactly once and predicted by a model fitted to every other family.

| Source model | Study-equal heldout RMSE | Improvement vs context | FAMILY-equal heldout RMSE | Improvement vs context |
|---|---:|---:|---:|---:|
| Original field context only | 0.29374 | baseline | 0.29789 | baseline |
| **Night/day bloom time** | **0.28146** | **+4.181%** | **0.26969** | **+9.466%** |
| Other 5 floral traits, NOT bloom time | 0.30172 | **−2.717%** | 0.30535 | **−2.505%** |
| Bloom time + other 5 floral traits | 0.29462 | **−0.298%** | 0.28749 | **+3.492%** |

A family has between **1 and 84** genuine published paired reproductive observations (median **3**). The largest family therefore accounts for **84/196 = 42.9%** of the original study site/year data. Only **8 of 28 botanical families contain both day-blooming and night-blooming plants** in the source. These are serious taxonomic opportunity and sampling limitations.

Of the 28 families, **18** have lower held-out squared error with bloom-only than field-context-only prediction. Fixed random-seed 2,500-replicate FAMILY bootstrap of the family-equal RMSE gain gives a **descriptive 95% interval [−3.43%,+22.20%]**. The analogous combined bloom+other-morphology range is [−12.84%,+18.85%], and morphology alone [−15.24%,+4.85%]. These intervals are over selected published families, not a random sample of the entire plant kingdom. They DO NOT establish universally transferable night-bloom performance.

## 48 leave-one-whole-genus-out predictions

Across the same 196 original observations, with source FAMILY equal total weighting: context-only RMSE **0.29811**; bloom-only **0.27234 (+8.642%)**; other traits-only **0.30475 (−2.228%)**; bloom+other traits **0.28948 (+2.894%)**. Under original study-equal weighting, bloom-only improves **2.954%**, other traits-only **−2.915%** and combined **−1.425%**. These predictions hold entire genera out, but may include the same genus' study context if other genera were tested by that paper; no randomized plant phenotype manipulation.

## Five-fold held-out published studies, ALSO purge every held-out plant species

Simply holding out papers can train on the same species from other papers. This stricter evaluation drops **13–56 additional training rows per fold** with plant species present in the held-out original studies. All 5 folds retained adequate independent training study support, with predictions on identical source test cases:

| Fixed original model | Study+species-purged RMSE | Improvement vs field-context baseline |
|---|---:|---:|
| Field context only | 0.29666 | baseline |
| **Bloom time only** | **0.28065** | **+5.396%** |
| Other five floral traits only | 0.31070 | −4.732% |
| Bloom time plus five other traits | 0.29421 | +0.826% |

Thus blooming time contributes some out-of-study/novel-species functional information; generic syndrome morphology does not provide additional stable held-out predictive benefit.

## Original night-bloom functional association if each botanical FAMILY is fully omitted

All **28 of 28** leave-one-family-out source study-cluster WLS fitted successfully. The night-bloom versus day-bloom adjusted paired reproductive share coefficient ranged **+0.16494 to +0.25409 fractions** (+16.494 to +25.409 pp), with all **28 nominal original-study-cluster CIs completely above zero**. This is stable point sign but those 28 estimates are highly correlated refits of a single database, NOT 28 independent experiments. Only **4 original studies** compared day- and night-blooming plants within one paper and only 8 botanical families compared both classes. Taxonomy and paper variation cannot be fully separated.

## Previous five-fold family transfer and important reconciliation

The original PR 54 five-fold held-out FAMILY analysis gave blooming time study-weighted RMSE gain only **+0.146%** versus field-context baseline, morphology alone **−6.59%**, and bloom+other traits **−5.75%**. The new **28-fold leave-one-complete-family-out** study-weighted bloom-only gain **+4.18%** (family-weighted +9.47%), but its original FAMILY bootstrap 95% range still includes 0. The difference matters: training-cohort composition changes across folds and the original study/family distribution is highly imbalanced. We do not portray any one fold strategy as proven universal predictive performance.

## Scientific conclusion / what has NOT been established

**Supported in this selected real experimental published source:** Night-blooming species are associated with greater actual night-only (relative to day-only) fruit or seed production within documented day/night pollinator-exclusion trials; the association persists whole-family deletion, while bloom timing alone offers SOME transfer to unseen original source studies/species/genera/families on certain fixed protocols and weightings.

**Not supported:** a universally transferable joint color+shape+scent+nectar+length floral syndrome, particularly beyond known families. Published family morphology-only trait models worsen held-out predictions. This does not mean odour or UV/color cannot be essential at a particular plant species' true pollinators.

**Not demonstrated:** actual local moth-specific pollen deposition, an ALAN manipulation causing effective pollinator decline, an island-only plant functional shift, autonomous selfing compensation, or evolutionary simplification. Original field light contrast is primarily between published studies. Study selection and trait measurement are external and observational. Missingness, family imbalance and family bootstrapped uncertainty are material limitations.

## Provenance and executable source

Source pipeline: scripts/analysis/alan_real_DvN_flower_function_family_LOFO_purgedCV_20261113.py; tests/test_alan_real_DvN_flower_function_family_LOFO_purgedCV_20261113.py; workflow .github/workflows/alan-real301-floral-function-family-transfer-20261113.yml. Independent original 301 paired reproduction and 30 signed annual PCNL raster sources from run 43380089481, and prior same-source complete floral trait PR 54, run 43401490182. Initial source result run 43402398580 SUCCESS; FINAL training-only study-weight recomputation run 43403060451 SUCCESS. Bootstrap/metrics above use the final validated rerun. Only aggregate model results and count summaries are published; raw author site/species/source trait-by-location data stay in original controlled parent data. PR 55 is Draft, no previous ALAN results or main have been changed.
