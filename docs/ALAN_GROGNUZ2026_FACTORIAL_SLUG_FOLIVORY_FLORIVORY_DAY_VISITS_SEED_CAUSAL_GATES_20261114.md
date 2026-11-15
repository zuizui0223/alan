# ALAN: does light-mediated herbivory truly compensate reduced pollination? Two distinct damage pathways

**2026-11-14, independently audited PUBLIC PREPRINT evidence, not source-row reanalysis.** Source Grognuz & Knop (2026), *Artificial light at night restructures daytime plant-pollinator interactions through direct and herbivore-mediated pathways*, [public author-uploaded original March 2026 preprint](https://www.researchgate.net/publication/401424487_Artificial_light_at_night_restructures_daytime_plant-pollinator_interactions_through_direct_and_herbivore-mediated_pathways), DOI [10.13140/RG.2.2.36822.10563](https://doi.org/10.13140/RG.2.2.36822.10563), related University of Zurich thesis [10.5167/uzh-434379](https://doi.org/10.5167/uzh-434379). Original source authors' [2026-04-10 Mendeley Data 10.17632/drhzznh4ps.1](https://data.mendeley.com/datasets/drhzznh4ps/1) lists the underlying study files. **Mendeley original sample-level data were not acquired or reanalysed here; the API had previously returned 401.** Do NOT mislabel preprint as peer-reviewed published or re-computed original source regression.

**Evidence-path CI/sign audit SUCCESS:** [GitHub Actions #43574842738](https://github.com/zuizui0223/alan/actions/runs/43574842738), [aggregate source-number receipt artifact #14106076319](https://github.com/zuizui0223/alan/actions/runs/43574842738/artifacts/14106076319).

## 2026-11-15 genuinely executed new source-access audit: original data still inaccessible without authorized access

The published, CC BY 4.0 [Mendeley source record version 1](https://data.mendeley.com/datasets/drhzznh4ps/1), DOI 10.17632/drhzznh4ps.1, is genuinely indexed and its public HTML landing page returns **HTTP 200** (118,538 bytes). But the actual individual experimental data must not be claimed downloaded. A new source-probe workflow tested the published Mendeley public data/snapshot/file-list APIs and documented **401** from each, the supposed public zip route returned **404**, and the source page HTML had **no directly recoverable UUID entries**. **Zero original raw research files downloaded; zero new raw-seed regressions fitted**. Source access is a technical blocker, **not evidence that the study lacks data or that its ALAN effect is absent**. The already reported preprint light-to-seed beta −0.40 [−0.73,−0.08] remains **author-reported, not independent reanalysis**. [GitHub Actions #43682158500 SUCCESS (source access receipts only)](https://github.com/zuizui0223/alan/actions/runs/43682158500), aggregate-only artifact #14148369522, ZIP SHA256 63c75b918a78bef865ea375a92a9841e25c4290648280fe2ed3dd92a8df71721. Source code: scripts/analysis/alan_Grognuz2026_slug_2x2_public_raw_probe_20261115.py. No bypass of authentication was attempted and no treatment data fabricated. Alternative independent directly replayable reproduction evidence is the 3-field author source [Li 2026 DOI 10.64898/2026.07.10.737511](https://doi.org/10.64898/2026.07.10.737511), now separated into seed-positive hurdle and fruit/seed counts in [the 2026-11-15 original-source report](https://github.com/zuizui0223/alan/blob/analysis/alan-floral-timing-dependence-20261115/docs/ALAN_LI2026_REPRODUCTIVE_TWO_STAGE_HURDLE_20261115.md). The two experiments are independent, with different plants/treatments/outcomes, not replicates of one light×slug model.

## Real experiment reported by the authors

For two actual native Central European wild plants (*Trifolium pratense*, *Centaurea jacea*), original 2023 potted plants were exposed to **four 2×2 combinations** of 3 consecutive nights under lit/unlit artificial night light × nocturnal slug presence/absence. After exposures, source plants were combined into a common daytime field arena and naturally visiting insects were recorded in short standardized surveys; source paper reports **819 actual insect visits**, **28 survey rounds**, 20 original Trifolium pot blocks and 8 Centaurea blocks. The paper reports **276 herbivory-exposed potted plants**. Seeds were measured subsequently. This is ONE original factorial study with nested potted plant/block structure, **not 28 independent published experiments** or a multi-site tropical island survey.

## Actual source-reported effects (no re-fitting)

| Reported causal/statistical relationship | Original β or odds ratio | 95% interval | Original reported p | Direct empirical interpretation |
|---|---:|---|---:|---|
| Nightlight vs darkness: subsequent total daytime visits (herbivore-free context) | β −0.38 | [−0.60,−0.16] | .013 | Decreased visits |
| Nightlight: subsequent Hymenoptera visits | β −0.30 | [−0.52,−0.09] | .018 | Decreased bee/wasp visits |
| Nightlight: subsequent Diptera visits | β +0.06 | [−0.05,+0.25] | .49 | Unresolved |
| Nightlight: odds of slug-damaged leaves | OR **0.21** | [.045,.99] | .048 | 79% lower **damage odds**, NOT necessarily 79% fewer damaged leaves |
| Folivory (leaf damage) main association with visits | β −1.3 | [−2.6,+0.05] | .060 | Unresolved |
| Nightlight × folivory for visits | β +2.1 | [−0.27,+4.4] | .083 | Unresolved |
| Florivory (flower damage) association with visits | β **−3.9** | [−6.2,−1.7] | <.001 | Fewer visits with flower damage |
| Nightlight × florivory for visits | β **+4.6** | [+2.3,+7.0] | <.001 | Conditional visit association/effect modification, not guaranteed mediation |
| Nightlight and subsequent seeds | β **−0.40** | [−0.73,−0.08] | .014 | Lower original author-model seed outcome |
| Nightlight × herbivory for seeds | not numerically resolved in extracted model summary | — | Not significant as reported | No demonstrated compensation of seed outcome |

Original β units/model link scales across visit, interaction, and seed models are **not reconstructed from raw original R scripts**; never express β as percentage change or directly add these different-model coefficients. The original seed interaction is *not* proven zero; it was not statistically resolved.

## The previously underappreciated critical distinction: different kinds of herbivory

**LEAF DAMAGE PATH**
- The original lit treatment decreased slug leaf damage odds significantly (OR 0.21), but observed leaf damage → decreased daytime pollinator visitation has a reported CI that includes zero (β−1.3, CI −2.6 to +0.05). Likewise ALAN×leaf damage for daytime visits remains uncertain (β+2.1, CI −0.27 to +4.4).
- So ALAN → less folivory is documented, but the SAME leaf-damage mediated increase in daytime **effective pollen service** and mature seed set is **not established**.

**FLOWER DAMAGE PATH**
- Damage to floral structures has a strongly negative reported association with daytime visits (β−3.9), and the light×flower-damage interaction is positive (β+4.6).
- However, the authors report original **ALAN did not detectably change florivory incidence** (around 33.5% of 355 examined inflorescences damaged). The strong light interaction thus cannot be combined with the independently significant reduced **leaf** damage to claim a complete ALAN→less **flower** damage→higher visits pathway.
- Flower damage amount is a **post-treatment measurement** rather than a separately randomly assigned intervention in the 2×2 ALAN×slug-present factorial, so interpreting the reported light×observed-florivory coefficient as an isolated causal mediating mechanism requires care.

**SEED ENDPOINT**
- In the original preprint, light had a negative reported seed-set coefficient (β−0.40, p=.014), while herbivore main/interactions were not demonstrated to rescue seed set. Hence visit-level buffering and net seed reproductive compensation are distinct and the latter is **not established**.

## Ecological interpretation for the GLOBAL ALAN programme

The experiment supports the possibility of two opposing **trophic/time channels** after nighttime light (a direct negative daytime floral-attractiveness/visitor effect and changes in nocturnal slug damage) in the same original potted plant system. But **the two individually statistically supported damage edges do not form one complete causal mediation chain for seed formation**, because their damage variable types and endpoints differ. The causal light→floral VOC change is suggested, not directly verified, and no gecko/bat/spider has been demonstrated in these treatments to replace moth pollen transfer.

This is a distinct experimental pathway from:

1. Lu 2024 Tibet original 5-block day+night interactions: nocturnal Lepidoptera ↓, Diptera ↑, but no effective flower pollen or fruit on every source taxon.
2. Giavi 2016 Swiss six paired fields: daytime guild composition changes under light, but Swiss visits must be compared after original floral-cover normalization; our new source [Draft PR #69](https://github.com/zuizui0223/alan/pull/69) finds raw 5/6 sites positive for Diptera reduces to only 3/6 after plant cover adjustment.
3. Boom 2020 four-light Netherlands: effective pollen mutualist moth larvae consume seeds; fruit-set similar, light alters fruit infestation, estimated net seeds/flower model ANOVA p=.63.
4. Kays and Wenzell authentic day/night experiments: direct source animal groups and fruit outcomes, **without manipulated light**; cannot be combined as light-causal strata.
5. 562-island global accepted-photo flora: enormous photographed *opportunity* but only 52 directly externally night/day reproductive-characterized source species; no same-island ALAN×predator×phenology matching.

These are distinct outcomes. Do not produce a synthetic “four-study combined pollen service percent change” or assign established island/global generality based on 2 native Swiss potted species.

### Testable macro-ecological design changes

For future original light experiments and comparative observational island/mainland programmes, use separate standardized variables per study site×year:

- ALAN lux/spectrum and source exposure, matched lit-dark original field block; observational VIIRS/PCNL kept separate.
- Observed **leaf damage** vs observed **flower damage** and their nocturnal animal identities, alongside true predator attack and **predator pollen service** as separate paths.
- Recorded daytime/nighttime visits, direct stigmatic pollen and mature net seeds, with flower abundance and flowering time source profiles.
- Realized pollination dependence and autonomous selfing; classify existing island moth shift from independent ancestral community support rather than color/flower tube proxy.
- Estimate mediation only when ALAN, THE SAME damage mediator, effective pollen delivery and plant fitness are all measured at the correct plant/site/time grain with adequate independent original studies.

## Reproducibility and peer-review state

Machine-readable preprint reported-effect ledger: config/alan_public_Grognuz2026_herbivory_daypollination_seed_preprint_20261114.json; source-value and CI sign validator: scripts/analysis/alan_Grognuz2026_preprint_ALAN_slug_visitor_seed_evidence_gates_20261114.py and associated synthetic source tests; GH Actions workflow alan-Grognuz2026-slug-visit-seed-preprint-path-gates-20261114.yml. Run #43574842738 SUCCESS with aggregate source annotation, not raw table effects. Original manuscript read and full text source-linked, **no author raw data retrieved and no new effect-size model**. Main and all prior PRs unchanged; this Draft stacks on PR #69.
