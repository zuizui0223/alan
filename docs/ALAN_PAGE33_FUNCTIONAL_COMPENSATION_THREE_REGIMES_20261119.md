# ALAN mechanism Result B: three regimes of static functional compensation

2026-11-19 | Source-verified Page original 623 taxon rows; **new reanalysis / not an ALAN treatment** | Draft only

## What is genuinely new here

In independent empirical *single-visit conspecific-pollen-delivery* networks, how often is numerical pollinator replacement sufficient in principle to restore pollen delivery when the currently most frequent visitor taxon is lost?

This evaluates a concrete, falsifiable premise of ALAN plant-pollinator network interpretations: "new visitors replace lost visits, so pollen supply is restored." The *numerical sufficiency* of replacement depends on source-measured per-visit effectiveness; it cannot be read off visit counts alone.

The source data are Page et al. (2021), DOI 10.1002/ajb2.1764, original Figshare 10.6084/m9.figshare.16160265.v1 (623 taxon-level lines, source SHA256 e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a). Exact original study × plant × period × source metric and visitation units yield 72 measured source networks from 33 original studies and 62 plant names, **not 72 independent lighting installations**. The original literature was assembled for a different meta-analytic question (honeybee single visit effectiveness), not as a representative sample of ALAN-exposed wetlands.

## Frozen counterfactual and exact arithmetic

For each original network, let V_j be taxon j's source-reported visit rate and q_j its source-reported mean single-visit conspecific pollen deposition. Its initial pollen-delivery *proxy* is

P = sum_j (V_j q_j).

Select j* = argmax_j V_j (the most frequent recorded visitor, with deterministic source order for ties). Remove all visits of j*, **without changing the observed q**. Reassign precisely V_j* new visits to remaining taxa in one of two hypothetical ways:

1. **Proportional replacement**: distribute these visits in proportion to each remaining taxon's original visit rate, so
   P_prop = sum_(k != j*) V_k q_k + V_j* [sum_(k != j*) V_k q_k / sum_(k != j*) V_k].
   Define the required count multiplier kappa = q_j* / qbar_remaining. If kappa > 1, equal-count proportional replacement is insufficient.
2. **Best-case replacement bound**: allocate *every* lost visit to the most effective taxon among remaining, so
   P_best = sum_(k != j*) V_k q_k + V_j* max_(k != j*) q_k.
   If P_best < P, even optimal count-conserving reallocation fails **under this static network assumption**.

The scenario does **not** measure pollinator removal, actual visit redistribution, additional flowers, q changes after lighting, pollen competition, pollen viability, or seeds.

## Newly executed empirical counterfactual

[Original-source Github Actions SUCCESS](https://github.com/zuizui0223/alan/actions/runs/44545544302); 8 synthetic arithmetic/source-gate tests and source-locked analysis succeeded. [Exact aggregate source-verified JSON](../results/all_species_forest/biotic_interactions/page_compensation_three_regimes_verified_20261119.json).

| Counterfactual response class | Networks | Proportion | Original study-cluster descriptive 95% interval |
|---|---:|---:|---|
| Even the *best* remaining visitor cannot compensate with the same number of visits | **21/72** | **29.2%** | [17.7%,43.2%] |
| Proportional replacement fails, but best-case targeted rewiring could compensate | **18/72** | **25.0%** | [13.0%,36.0%] |
| Original proportions of remaining visitors are sufficient for equal-count compensation | **33/72** | **45.8%** | complementary partition, no independently estimated interval |
| Proportional replacement insufficient, summed | **39/72** | **54.2%** | [39.2%,69.2%] |

Across all 72, the mean positive-only pollen deficit under *proportional* full-visit replacement is 14.56% [9.76,19.74]; under the hypothetical *optimal* full-visit replacement 7.50% [4.22,13.07]. These are averages of **max(0,1−P_scenario/P_original)**, not net percent decreases or observations from field lighting experiments. The median *finite* required visit multiplier under proportional replacement was 1.077, and 15/72 (20.8%) required over twice the lost visits.

### Sensitivity

- At least five original visitor taxa: 43 networks /16 source studies (42 after correcting the verified butterfly taxon). Proportional replacement fails in **26/43**, optimal fails **10/43**.
- All positive original single-visit q: 58 networks /30 studies in the uncorrected original. Proportional fails **33/58**, optimal fails **19/58**.
- Exclude the original mislabelled butterfly *Vanessa carye* from the moth category: the main full-sample three-state counts **21,18,33** are unchanged. Do not pretend all original insect taxa are independently taxonomically reverified.
- **Moth most visited**: only 5 networks from 2 source studies; none required more replacement visits in the frozen scenario. This cannot prove that moth replacement is generally easy. It is a crucial counterexample to a universal "moth loss always functionally irreplaceable" framing.
- **Fly most visited**: 6 networks from 4 source studies, including 1 optimal-replacement failure; cannot estimate an independent fly-specific meta effect.

Study-cluster bootstraps are descriptive uncertainty, not significance from 72 independent experiments. Missing source q sampling SE and removal-induced q changes are outside the intervals.

## Can we connect functional compensation to the observed ALAN wetness bias?

**Not with the present original studies.** An exact-name audit against independent European Tichy EIV MOIST (official source SHA256 ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba) yielded:

| Source | EIV score coverage | Wet MOIST >7 | Within-study wet+dry support |
|---|---:|---:|---:|
| Page 72 networks /62 plants | **28 networks /25 plants** | **3 networks /2 plants** | 3 original studies both wet and dry |
| Lu 2026 field experiment's 11 focal species (5 actual light pairs) | **1 plant: Geum aleppicum (MOIST 6.0)** | **0** | unavailable |

Other missing traits are NOT scored dry; European EIV is not automatically transportable to the Qinghai/Tibetan alpine flora. No coefficient on *MOIST × ALAN* or *MOIST × functional threshold* is estimated here.

Independent wet-habitat macro finding (+0.534 adjusted MOIST difference, 158 matched candidate species) remains the main descriptive result in the original manuscript; nothing in this reanalysis confirms or falsifies its specific causal mechanism. The 2026 Lu five lighting pairs show altered day-night guild visits, but their actual plants do not match Page's five verified within-plant moth/fly SVE cases. Do not transfer q across species to build a fictitious observational mediation model.

## Article Results structure if this remains in the same manuscript

- Result 1: original light-associated moist-plant recording-bias finding and regional sensitivity (no change).
- Result 2: *experimental stage evidence*, source-authentic day/night visitor guild reallocation (Lu 5 real paired lighting units; already published, new reanalysis cannot masquerade as newly discovered average causal effect).
- **Result 3: newly inferred bounds on functional compensation**, the 3-state 21/18/33 distribution above. Label it "independent static pollination-network analysis", not "ALAN-caused reproductive failure."
- Result 4: wetness-by-mechanism trait overlap **NOT ESTIMABLE** in Page and Lu, in supplementary source audit.

Whether Results 2–3 fit the main manuscript ultimately depends on a future *independent* ecological bridge. To distinguish direct plant light responses, pollinator substitution, and wetland drainage, the next original experiment needs matched same-plant lit/dark ΔV by taxon, measured q or compatible pollen outcomes, wetness traits for both wet/drier plants, and plant reproductive outcomes, with independent replicated lighting installations. Preselect source-year, site and taxon overlaps before modelling.

## Reproduction

- [Counterfactual arithmetic, source and original-study bootstrap](../scripts/analysis/alan_page2021_replacement_threshold_20261119.py)
- [Source-locked EIV plant trait and Lu 11-plant exact-name support](../scripts/analysis/alan_page2021_moisture_function_overlap_20261119.py)
- [Unit tests for proportional/best-case bounds](../tests/test_alan_page2021_replacement_threshold_20261119.py)
- [Main CI workflow](../.github/workflows/alan-page-functional-replacement-20261119.yml)
- [Complete 72-network per-source artifact](https://github.com/zuizui0223/alan/actions/runs/44545544302)
- [Summary including all original-study sensitivity bootstrap intervals](../results/all_species_forest/biotic_interactions/page_compensation_three_regimes_verified_20261119.json)

**Release gate:** keep draft. A generalized ALAN-to-wetness causal claim is STOP until original source matches or an independent experiment supports it. No main merge.
