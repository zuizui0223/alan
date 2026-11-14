# ALAN 2026: original field × species × flowering round matched light-arrangement audit

## Scientific aim
Can the Li et al. (2026) **positive response of seed-bearing plants under 26-bulb linear arrays versus unlit plots** be attributed specifically to light *geometry*, rather than field, plant species, flowering schedule, initial flower-bud number, light footprint and intensity?

Original preprint: [Li, Adeniji, Meah & Clements (2026), bioRxiv 10.64898/2026.07.10.737511](https://doi.org/10.64898/2026.07.10.737511) (not established peer-reviewed). The actual archival original `Dongboli/experimental-data/Artificial lights at night/seedset.csv` is fetched directly by SHA-locked code in PR #73 and analyzed in this new independent Draft. We do not redistribute original per-plant records.

## Raw-source denominator and design audit
- Original author preprint Methods describes **90 potential plant deployment slots**: 54 *Matthiola longipetala* over three consecutive two-week rounds; 18 each *Silene nutans* and *Zaluzianskya capensis* in a single round.
- The **published original author CSV actually contains only 79 source plant rows**, precisely **44 ML, 18 SN, 17 ZC**; 78 have nonmissing seed-presence outcomes. This discrepancy does **not** establish 11 plant failures or data misconduct. It requires archive/protocol reconciliation before assuming all planned plants were measured.
- Three genuinely independent experimental hayfields, with spatial array / point / unlit plots per field; positions were changed halfway through season. There are **NOT 79 independent lighting installations**.
- Crucially, the article's original Methods specify a **26-lamp 25 m array** versus a **4-lamp total endpoint point treatment** (2 per endpoint), with no effort to standardize total lamp brightness or illuminated footprint. A geometry-only estimand is therefore **not identified**.
- The source `pollinated` flag is operationally **at least one observed seed**, not stigma pollen delivery, and the actual nocturnal visitors were not recorded.

## Observed original outcomes already source-replayed in PR #73

| Actual arm | Seed-positive original plants | Source seed counts | Source total fruit counts |
|---|---:|---:|---:|
| Dark C | 8/26 (30.77%) | 672 | 32 |
| Point P | 12/25 (48.00%); 1 outcome missing | 967 | 33 |
| Linear L | 20/27 (74.07%) | 3,729 | 149 |

Within the **three independent fields**, linear-minus-dark **risk differences** were +0.3333, +0.5051 and +0.4556, all three positive. But crucial *direct linear-minus-point* differences were **−0.1667, +0.3333, +0.6143** (only two of three fields positive), matching the original preprint's explicit no-significant L-versus-P test.

More stringent original species×field direct L-versus-P comparison finds **4 positive / 4 tied / 1 negative** among 9 supported cells; no claims of nine independent field replicates. In contrast, the L-versus-dark plant-level fraction is positive across all three fields and six of eight source-comparable species×field strata (two tied).

## New source-diagnostic sensitivity: original initial flower-bud opportunity
The original author paper states that starting flower-bud number did not differ significantly between treatments (source paper p>0.05); we do **not** contradict that inferential claim. Direct raw CSV shows the descriptive mean initial buds per observed outcome-known plant approximately:
- **L 17.89** initial buds/plant;
- **C 14.19** initial buds/plant;
- **P 16.32** initial buds/plant.

For a noncausal availability-denominator sensitivity (NOT actual per-open-flower pollen efficacy), divide summed seeds/fruits within each treatment by summed **initial flower buds**:
- Linear seeds per initial bud **3729/483 = 7.72** versus dark **672/369 = 1.82** and point **967/408 = 2.37**.
- Linear fruits per initial bud **149/483 = 0.3085**, versus dark **32/369 = 0.0867** and point **33/408 = 0.0809**.

This baseline-opportunity normalization still preserves a positive observed source linear-control contrast, but **is NOT a matched causal adjustment**, because initial buds are not ovules or the complete flowering resource availability. It cannot identify nocturnal insect function or linear illumination geometry.

## Actual newly executed matched analysis (author-source replay)
[Initial signed original Git SHA replay #43582105435 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43582105435); all four synthetic guard tests and original source matched-source reanalysis passed. Aggregate-only artifact **#14109225001**, SHA256 `22f416d2f80db0284f56365f24a8b85cce37ac2284db4c41f590658469546a2d`.

| Original matched site-support | L-C | L-P | P-C |
|---|---:|---:|---:|
| Species × hayfield support (not independent sites) | 8 | 9 | 8 |
| Positive / tie / negative | **6 / 2 / 0** | **4 / 4 / 1** | **4 / 3 / 1** |
| Equal-support within-stratum seed-presence difference | **+0.4451** | **+0.2722** | **+0.2014** |
| Species × hayfield × round support | 12 | 13 | 12 |
| Positive / tie / negative in round-matched strata | **7 / 5 / 0** | **7 / 4 / 2** | **3 / 6 / 3** |
| Equal-support within-stratum risk difference | **+0.4167** | **+0.3462** | **0.0000** |

These 8–13 cells are **not 8–13 independent experimental hayfields**. They represent the same three original spatial hayfields and deliberately do not report misleading cluster p-values.

Sensitivity to the flowering round: L versus C is **+0.4792** at round 1 (all three taxa), **+0.8000** at round 2 (only ML, n=3 L and n=5 C), and **+0.0333** at round 3 (only ML, n=6 L and n=5 C). The species composition and round are confounded; there is no valid general season-dependent ALAN coefficient. L-P is **+0.1042, +1.0000, +0.3333** by rounds, with extremely small round-2 samples. P-C changes sign across the rounds: **+0.3750, −0.2000, −0.3000**.

Leave out each entire **plant species** (not independent geographic replication): source pooled array-minus-dark seed-presence differences remain positive **+0.5333** (without ML), **+0.4143** (without SN), and **+0.3961** (without ZC). The positive array-control direction is not explained by a single plant species in these descriptive subset checks.

### Worst-case bounds for the 90 potential versus 79 archived plant observations
Assume for sensitivity ONLY that the 11 unavailable potential source rows represent unknown outcome-bearing individuals. Allocate any hypothetical missing plants in the **most unfavorable possible way**: extra **failures** to the first arm and extra **successes** to the comparison arm (remaining treated arms can be ignored). No real missing trial treatment assignments are asserted. For each number M from 0 upward, exhaust all integer allocations x failures to A and M−x successes to B; minimize `p_A-p_B`.

| Comparison | Observed source risk difference | Worst possible allocation of all 11 unrecorded plants | First possible algebraic reversal |
|---|---:|---:|---:|
| L−C | **+0.43305** | **+0.18671** (6 L failures, 5 C successes) | **23** such hypothetical missing observations |
| L−P | **+0.26074** | **+0.03571** (8 L failures, 3 P successes) | **14** hypothetical |
| P−C | **+0.17231** | **−0.04202** (3 P failures, 8 C successes) | **9** hypothetical |

This calculation is a **mathematical worst-case pooled-direction bound**, not a scientifically established missingness model, proof that potential deployments occurred, or field-/plant-stratified causal inference. In particular, the per-field signs may be changed by allocating unknown observations to specific small fields; no claimed field-robustness to the 11 absent potential rows.

## Ecological interpretation
The strongest source-replicated claim is a **conditional POSITIVE reproductive response to one compound linear-array treatment compared with an unlit control**, robust in direction across three fields, and preserved in observed source seed/fruit/bud ratios.

Claims **NOT** supported:
- linear light geometry independent of total photons/illuminated area;
- higher effective moth pollen delivery;
- observed moth corridor/habitat connectivity or predator shifts;
- specific ALAN-induced loss or gain on islands;
- a significant array-versus-point effect across independent fields.

Other original ecological anchors remain separated: Knop et al. (2017) 62% lower nocturnal visitation and 13% lower focal fruit set; Lu et al. (2026, Ecol Evol) 57.4% nocturnal Lepidoptera visit reduction with no measured seed in that network experiment; Lu et al. (2026, J Plant Ecol DOI 10.1093/jpe/rtag160) independently reports 56% lower seed set for nocturnally dependent self-incompatible *Pedicularis kansuensis* under artificial light, but little change for autonomously selfing *Euphrasia regelii* — published-author reports, **NOT new raw-data reanalysis**. Collectively these suggest testing `ALAN exposure regime × independently measured nocturnal reproductive dependence × autonomous selfing`, not a single global directional light slope.

## Reproducibility
Real SHA-pinned original file, original species/field/round matching with explicit missing-arm NAs, source initial-bud normalization, leave-one-species-out and same-round diagnostic are in `scripts/analysis/alan_Li2026_original_matched_strata_bud_audit_20261114.py`, four mathematical tests in `tests/test_alan_Li2026_original_matched_strata_bud_audit_20261114.py` and native GitHub Actions `.github/workflows/alan-Li2026-matched-field-species-bud-20261114.yml`. Output only aggregate JSON; no raw identifiers/coords deposited.
