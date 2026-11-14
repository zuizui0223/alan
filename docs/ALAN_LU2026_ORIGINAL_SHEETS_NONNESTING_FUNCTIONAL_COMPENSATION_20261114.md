# ALAN: source-original Lu (2026) cross-sheet reconciliation and functional-compensation sensitivity

## Purpose and evidence level
This independent audit reuses the verified original published Lu et al. (2026), *Ecology and Evolution*, DOI 10.1002/ece3.74368, supplemental `Date.xlsx` (SHA256 `ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f`). It re-downloads the authentic workbook during CI and runs the existing frozen source parsing from PR #65, without manufacturing extra replicate treatments. Exactly **5 source field blocks**, **400 survey-level visit rows**, and **1,038 separate focal plant-visitation rows**. Neither worksheet measures pollen deposition/seed set.

## Newly flagged source-comparability failure
The 11 focal plant species visitation worksheet has:
- Dark nocturnal Diptera: 5 moth-majority plants **87** plus 6 others **19** = **106** recorded fly visits, but the nominally all-survey `visitation` worksheet has only **54**.
- Lit nocturnal Diptera: **133 + 48 = 181** focal fly visits, but all-survey has **151**.
- Difference focal minus all-survey is **+52 dark** and **+30 lit**.

Therefore the two published author worksheets **cannot be treated as a literal same-unit subset/partition**, regardless of any rounding or composition model. Other insect orders need not show the same discrepancy. This is **not proof that the authors' experiment or original observations are incorrect**: worksheets may use different counting rules, survey subsampling, per-plant multi-counting, or definitions. The required follow-up is a source-method/author reconciliation of event definitions; do NOT subtract focal rows from global survey totals or pool them into a meta sample.

## New conditional functional compensation calculation (SCENARIOS ONLY)
Within **each sheet independently**, let a source-recorded Lepidoptera visit have hypothetical pollen-delivery weight 1, and a Diptera visit weight `r`. Ignoring other orders, the illuminated/dark weighted-count ratio is:

`R(r)=(L_lit+r*D_lit)/(L_dark+r*D_dark)`.

The break-even point is `r*=(L_dark-L_lit)/(D_lit-D_dark)` when moth-labelled visits fall and fly visits rise.

| Source worksheet / stratum | Dark Lepidoptera→lit | Dark Diptera→lit | r* required to balance those TWO visitor orders |
|---|---:|---:|---:|
| All 400 original survey rows, nighttime | 838 → 357 | 54 → 151 | **4.959** |
| Distinct 1,038 plant-row worksheet: 5 nocturnal Lepidoptera-majority focal plants | 680 → 255 | 87 → 133 | **9.239** |
| Same plant-row worksheet: other 6 focal plants | 95 → 73 | 19 → 48 | **0.759** |

These are **algebraic requirements**, *not estimated biological fly effectiveness*. Values in the two worksheets are deliberately **not** merged. At r=1, the source-recorded Lepidoptera+Diptera weighted-count lit/dark ratio is 0.570 in the all-survey sheet, 0.506 for the five focal plant taxa, 1.061 for the other six. At r=5, the two-order ratio is about 1.004 (all survey) and 0.825 (five focal); all other guilds, per-flower availability, time-bin mismatch and true stigma pollen/seed remain unmeasured.

## Out-of-block robustness guard
In original PR #65's source-independent whole-experimental-block holdout, **only 4/5 blocks** have at least two positively documented focal plants in both groups and **only 2/4** support the pooled "Lepidoptera-majority loses more visits" ordering. Therefore the above subgroup threshold **cannot be generalized to an evolved moth floral syndrome or global islands**. The robust result remains treatment-linked nighttime moth-labelled visits down, fly visits up, total visits down in one genuine original experiment.

## Reproduction
Source original and frozen parser from PR #65 are imported; new script `scripts/analysis/alan_Lu2026_cross_sheet_functional_threshold_20261114.py`, tests `tests/test_alan_Lu2026_cross_sheet_functional_threshold_20261114.py`, and source-verified Actions workflow `.github/workflows/alan-Lu2026-crosssheet-functional-threshold-20261114.yml`.

Run workflow on this branch to obtain a signed aggregate-only JSON artifact. No raw published workbook or original plot/plant-survey data are copied into the repository.
