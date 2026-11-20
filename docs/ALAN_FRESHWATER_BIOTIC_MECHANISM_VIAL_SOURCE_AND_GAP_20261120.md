# ALAN wet-biotic mechanisms: original data comparative feeding Result and meta-source boundaries (2026-11-20)

## Evidence classification

A result may appear in the ALAN paper's Results only with an explicitly stated estimand, original observational/experimental units, source-authentic uncertainty, and a new analysis question. Previously published effects are background evidence, not our new discovery. The macro overrepresentation of habitat-wet plants among long-term recording-decline candidates remains the first Results observation. **No analysis below establishes that aquatic consumer responses mediate this MOIST affinity bias.**

## Original three-consumer experimental data now accessible and reanalyzed

Campos et al. (2026), Frontiers Ecol Evol DOI [10.3389/fevo.2026.1830790](https://doi.org/10.3389/fevo.2026.1830790). Official original data deposited [Zenodo 19022243](https://doi.org/10.5281/zenodo.19022243) as \`data_Campos et al 2026.xlsx\`; SHA-256 \`b5c36d9b4c184a3f563388d56a9072fcfe9fa8c988600e9979054402083a7d6d\`. Six original sheets; three used here are feeding rates of freshwater consumers of conditioned (dead) alder leaves. Experimental unit in original methods is the **vial**, not each of its individual organisms. Our source-verified original-data rerun and tests [GitHub Actions #44660510047](https://github.com/zuizui0223/alan/actions/runs/44660510047) passed.

### New original-source ratio-of-means reanalysis

Within each source species, divide the mean vial feeding rate at 1 or 10 lux by the mean rate at 0 lux. Percentile intervals resample original vial rows *within light level* 4,999 times; they **do not** include between-light-installation, source-study, or cross-species model uncertainty.

| Consumer (source vial normalization) | Vial 0/1/10 lux with feeding value | 1 lux / 0 lux | 95% vial-bootstrap interval | 10 lux / 0 lux | 95% interval |
|---|---|---:|---|---:|---|
| **Potamopyrgus antipodarum** (per surviving animal; 7 days) | 10 / 10 / 10 | **1.334** | **[1.078,1.624]** | 1.233 | [0.994,1.521] |
| **Sericostoma vittatum** (per mg wet animal mass; 7 days) | 20 / 19 / 18 | 1.092 | [0.865,1.350] | 1.003 | [0.729,1.292] |
| **Chironomus riparius** (per surviving animal; 14 days) | 15 / 14 / 14 | **0.502** | **[0.358,0.673]** | **0.424** | **[0.291,0.591]** |

A new *between-taxon contrast in the ratio of response ratios* at 1 lux (P. antipodarum relative to C. riparius) is **2.658**, vial-bootstrap 95% [1.824,3.976]. At 10 lux it is **2.911** [1.951,4.518]. These are not independent-research estimates of diet/habitat moderation: experimental duration and source species differ (7 vs 14 days), and the original denominator differs for S. vittatum. None of these results tests actual grazing on live wet plants or their reproduction. Species-level directions agree with the authors' reported significance pattern: snail increased feeding only at 1 lux, chironomid decreased at both, caddisfly no clear change.

Original report published 20 vials/treatment for S. vittatum, 10 for P. antipodarum, 15 for C. riparius. Missing vial outcomes (S 0/1/2, C 0/1/1 across 0/1/10 lux) were not invented as zeros. For feeding, the original article reports no ROUT technical-outlier exclusions; its reported technical-outlier removals concern biochemical responses (S. vittatum GST/tGSH, P. antipodarum lipid peroxidation).

### Interpretation and why not a main wet-plant mechanism

These three consumers feed on **previously conditioned dead Alnus leaves** under laboratory illumination, not a living plant's leaves and not original wet-habitat decline candidate species. Therefore these empirical figures establish taxon- and dose-dependent ALAN alteration of **aquatic detrital consumption**, not ALAN → living wet-plant damage → plant loss. These three response traits are from one original laboratory study, with 0/1/10 lux level, not three independent experiments in a random-effects meta-analysis. Our confidence intervals assume the original vial is the independent replicate, following the published analysis; unreported shared fixture structure could narrow the effective design.

**New ecological question supported:** Does exposure target (plant vs consumer) and the consumer feeding mode predict direction of trophic effects? A future multi-experiment meta-regression needs independent ALAN manipulations and a harmonized feeding endpoint, explicitly separating live macrophyte grazing, conditioned litter, flower visitation, predator interception and reproductive outcomes.

## Separate live-macrophyte original experimental findings — existing published results only

- [Mondy et al., 2021](https://doi.org/10.1016/j.aquabot.2021.103447): *Ceratophyllum demersum* illuminated during growth, snails tested without ALAN. Original paper found 1.6× snail consumption of light-exposed leaves, isolating a plant-mediated palatability route. This is **not** ALAN-specific nocturnal attack by the snail.
- [Xue et al., Ecology 2026](https://doi.org/10.1002/ecy.70458): manipulated light × invasive *Pomacea canaliculata* grazing in live native aquatic macrophyte mesocosms. Reported native community biomass suppression 48.0% without ALAN versus 87.2% under ALAN, while *Myriophyllum aquaticum* was avoided by the invasive snail. This is **published authors' result**, not a newly estimated independent effect or a paired plant-species MOIST contrast. Light treatment has only 3 independently allocated experimental zones per arm in the monoculture experiment, so treating every shoot as independent installation would be invalid.
- [Grognuz et al. 2025](https://doi.org/10.1016/j.isci.2025.112770): terrestrial slugs decrease feeding/activity under lighting. Different system, plant quality and consumer exposure; cannot average its leaf damage as the same effect as live aquatic plant biomass or litter decomposition.
- [Grognuz & Knop 2026 preprint](https://www.researchgate.net/publication/404486765): field ALAN daytime floral visitors, slug herbivory and plant outcomes, but may share experimental site/installation lineage with prior Swiss investigations; do not count as independent replication without exact light-site audit.

## Existing herbivory-pollination 171-study meta: full raw effects not present in public supplementary archive

[Haas-Desmarais et al., 2025 online / Annals of Botany 2026](https://doi.org/10.1093/aob/mcaf258) reports 171 original studies and 1,348 observations concerning herbivore effects on flower traits, pollinator attraction and reproduction, without ALAN intervention. We fetched the original EuroPMC supplementary ZIP (SHA256 \`530c78af0f1b1a8729b402f8b1ef457373ea8457c8b26c75f4011cf678940b07\`), inspected all nested files and source CRCs; [verified Action #44659515692](https://github.com/zuizui0223/alan/actions/runs/44659515692). It contains **Appendix S1.docx** (plots/sensitivity, 0 structured Word tables) and **Appendix S2.docx** (8-row summary moderator test table). **There is no original 1,348-effect data table** in this public supplementary archive, so a MOIST moderator cannot be computed from it without a separate original data deposit. [Original source gate JSON](../results/all_species_forest/biotic_interactions/haas_desmarais171_original_supplements_structure_20261120.json).

## Current manuscript-ready inference boundary

1. Retain the independently documented main macro anomaly: wet-habitat affinity is higher in species selected as recording decliners in increasingly lit sites (+0.534 adjusted MOIST units).
2. Real ALAN experiment: nocturnal moth visitation loss and fly visitor reallocation (Lu 2026, five paired lighting blocks), but source-block exact overlap test **15 observed vs 16.05 expected**, two-sided p=.402, therefore no preferential plant-level fly targeting.
3. Original Page observational functional networks: visitation ≠ effective pollen. Static best-case count replacement failed in 21/72 networks, but this failure was not significantly greater than a within-network random q assignment null. **Not ALAN and not matched to wet-plant taxa**.
4. Original wet-plant *biological mechanism* remains unidentified: aquatic plant-grazer interactions are mechanistic candidates with source experiments, but there is no same-plant, same-context bridge from a wet MOIST value to ALAN-dependent live plant outcome across enough independent species.

### Reproducibility
- [Original Zenodo 6-sheet schema](../scripts/acquisition/alan_freshwater_zenodo_original_workbook_schema_20261120.py)
- [Vial mean ratio and cluster bootstrap source analysis](../scripts/analysis/alan_freshwater2026_three_consumer_vial_feeding.py)
- [Vial synthetic tests](../tests/test_alan_freshwater_three_consumer_vial_20261120.py)
- [Source-verified output JSON](../results/all_species_forest/biotic_interactions/campos2026_three_freshwater_consumer_feeding_20261120.json)
- [Original source 171-study appendix audit](../scripts/acquisition/alan_herbivory171_original_supplement_20261120.py)

**DO NOT merge into \`main\` as mechanistic explanation of wet-plant recording bias. Preserve as independent supplementary or separate Results with explicit non-causal ecological scope until a genuine live-plant cross-species moisture×ALAN interaction dataset is obtained.**
