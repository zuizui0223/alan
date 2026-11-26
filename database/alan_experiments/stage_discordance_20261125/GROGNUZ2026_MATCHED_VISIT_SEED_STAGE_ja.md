# Grognuz & Knop 2026: published pot-level visit and seed table, matched reproductive-stage audit (2026-11-26)

## Scope and new evidence

This is a new analysis of **printed Table S3 pot-level data**, not the original Mendeley data file or a new independent programme. The existing draft PR70 (analysis/alan-grognuz-herbivory-daypollination-seed-pathgate-20261114, config blob 7f52ebe92eb2e3078169c788db99dafba51734f7) only recorded reported coefficients and did not jointly extract pot rows. Source: Grognuz & Knop, March 2026 preprint DOI 10.13140/RG.2.2.36822.10563, https://www.researchgate.net/publication/401424487_Artificial_light_at_night_restructures_daytime_plant-pollinator_interactions_through_direct_and_herbivore-mediated_pathways, Tables S2/S3, Methods. Primary dataset 10.17632/drhzznh4ps.1, CC BY 4.0, landing page checked; network DNS fails in container, original archive not obtained. Publisher PDF linked by ResearchGate returned cache miss, so the transcription is from indexed text, not a visually independently verified PDF. Peer review status not confirmed.

448 pot records = 28 choice blocks × 16 pots, 4 slots per treatment. First 20 blocks T. pratense, last 8 C. jacea; 5 and 6 choice dates. 23 missing seed outcomes, not zeros. Each pot is a deployment record, not necessarily a unique individual. Authors reused plants; lamps and dark/illuminated positions recur at the same field. Slug presence was factorial, pots randomized to treatment, but distinct randomization of light at independent sites is not established. Three nights of exposure followed by a common 90-minute daytime choice trial, target flowers subsequently rebagged and mature developed seeds counted. Neither 448 pots nor 28 choice rounds nor 11 dates are independent lighting programmes. No new formally admitted meta-analysis effect.

Transcription file grognuz2026_tableS3_transcription_20261126.csv blob 433e67105583796668005904ba8383923e880218. Arm headings A-/H-, A-/H+, A+/H-, A+/H+ interpreted as D0,D1,L0,L1 in blocks of four columns. All eight reproduced seed means agree with printed rounded means within .005; some reported SDs do not agree exactly. All reproduced visitation means agree within .011 but total 822 versus Results' 819. A spaced entry '1 1 (0)' (240627-VG slot13) read as 11 visits/0 seeds, with a one-visit and whole-block-exclusion sensitivity. Do not silently change other values to force totals to match. Column-order/rounding/extraction discrepancies require the original CSV or a verified image table. The printed seed means are not enough to establish every cell is correct.

## Date-weighted contrasts, matched pots

Restrict both endpoints to the same observed-seed pots in each arm. Average four arms separately within each round, then equal-weight round summaries within dates. Average dates within species. Light menu includes both slug arms equally. RRs use ratios of arithmetic group means; delta scores on date summaries preserve stage covariance. Approximate t4/t5 intervals in JSON are **date variability references**, not true design-based confidence intervals or independently replicated lamp effects. All six species × light-menu results retained. Zero counts included; no pseudocount. No date-independence inference.

| species | menu | visits RR | seeds RR | seed/visit relative response |
|---|---|---:|---:|---:|
| T. pratense | both | .976 | .714 | .732 |
| T. pratense | slug absent | .895 | .805 | .900 |
| T. pratense | slug present | 1.075 | .638 | .593 |
| C. jacea | both | 1.088 | .822 | .756 |
| C. jacea | slug absent | .653 | .302 | .462 |
| C. jacea | slug present | 1.774 | 1.289 | .727 |

The C. jacea no-slug relative stage contrast is .462 (date-based nominal .268–.797); drop ambiguous round .451 (.250–.812), but missingness sensitivity and true light replication are limiting. T. pratense overall seed-positive fraction RR .812 (date-based .713–.925), whereas seed count among positive outputs ratio .880 (.519–1.490). T. pratense complete-round sensitivity (exclude every round with a missing seed): seed RR .591, seed/visit .637, seed-positive .792 (4 dates/13 rounds); C. jacea complete-round support fewer than 2 dates, STOP. Overall stage ratios are negative direction but date variability intervals include 1 for both species. No claim of universal seed-per-visit efficiency reduction, or light-caused mediation.

Slug presence modifies point responses in C. jacea: light × slug visit response ratio 2.716, seed response ratio 4.274 (date reference 1.260–14.501). T. pratense interaction ratios 1.201 and .792 respectively. This is a contrast of **assigned slug presence**, not the published post-treatment continuous damage covariate. The species differ in experiment dates and plant reuse; do not infer a confirmed species-specific herbivore-mediated rescue. The original reported adjusted seed model gives ALAN beta -.40 (95% -.73,-.08), with no clear damage interactions; model-adjusted and marginal means need not agree.

## Missing outcome and unobserved visit audit

Among 425 pot records with seeds measured, 171 had zero recorded visits; 50 of those had seeds >0. Total developed seeds 1479, of which 229 occurred on zero-recorded-visit pots. These are descriptive, not evidence of autonomous selfing: short observation window, bagging leaks, missed visits, or different flower histories may contribute. Thus dividing seed number by recorded visits as individual-level efficiency would discard real seed outcomes and be undefined at zero visits. Use the ratio-of-marginal-response-ratios only. Paper reports 3/39 bagged control flowers produced seed, all in C. jacea; controls distinct from target pots with zero recorded visits.

Missing seed sensitivity uses **scenario** 0 to maximum observed seed count within species (44 T, 19 C). This is not a biological upper bound, and the calculation uses pooled-pot means, not the date-weighted estimator above. Both slug arms pooled: T seed RR .620–1.083; C .298–2.277. No-slug T .651–.935; C .097–1.438. Thus the strongest C estimate can reverse under these missing-value scenarios; not confirmed. Preserve missing seed outcomes and no causal claim. Formal meta-analysis still cannot use date-level standard errors as lighting-unit variance.

## Direct access/reproduction source acquisition audit

Macgregor2019 primary article DOI 10.1002/ecs2.2550 (5 plots × 4 rotating experimental runs, 169 females, 417 flowers) rechecked. Original supplement ecs22550-sup-0001-AppendixS1.pdf returned inaccessible; White Rose original PDF screenshot cache miss. Full-night light increases pollination success, conditional seed count and mass no resolved difference; no light-by-caging interaction detected, low power. No numeric model contrast reconstructed from unreadable figures. Zenodo 1743917 code only per existing audit, network retry failed DNS. Wilson2021 DOI 10.1093/icb/icab010 (yucca moth, direct light gradient and skyglow separate) no original raw supplement acquisition. Li and Lu already audited; no new programme from duplicate records. No contacts sent to authors. The source Czarnecka 2025 prior archive not mounted, files search returned no match. No invented new source counts.

## Saved output and validation

JSON grognuz2026_matched_visit_seed_stage_20261126.json, blob aaaa7b92dc75c71231e0b7d35b1ef472c1282aab; standalone script scripts/analysis/alan_grognuz2026_matched_visit_seed_stage_20261126.mjs. 440 arithmetic identities rechecked from saved results, max absolute error 8.881784197001252e-16. Source data transcribed from published text, not independently double-extracted; standalone Node/CI not run. Main, prior 108-effect catalogue, exploratory synthesis results unchanged. Formal new k=0. This directly advances the joint visit-seed data extraction, but does not establish a second independent light × pollinator-access × reproduction experiment or a wet-plant macro mediator.
