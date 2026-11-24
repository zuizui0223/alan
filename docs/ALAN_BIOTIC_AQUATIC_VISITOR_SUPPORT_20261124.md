# Recorded aquatic-hoverfly support and alternative visitors (2026-11-24)

Exploratory follow-up, not a new main result. Frozen original source: `b780370:results/all_species_forest/wet_network/flower_sources.csv.gz` (blob 8d75a38af58c9fd63a1d1c99b3e69cf6e39f830e) and `species.csv` (blob 388a950643cda5769a41fbcecc8b4d4194448f7c). The calculation was executed against connector-retrieved original CSV bytes in JavaScript; the companion Node script is a local replay adaptation, not a separately executed repository test.

363 wet-cohort rows, 58 candidates; exclude all four rows with duplicated accepted names, leaving 359 unambiguous rows. 3864 flower-source rows, 29 distinct visitor labels flagged aquatic_hoverfly by the existing EuPollTrait join (on_water_plants or submerged_sediment); flags internally consistent. Do not interpret an unflagged visitor as terrestrial, or an unrecorded interaction as absent. No new taxonomy verification. Count distinct visitor names per plant, not number of records, visits or pollen grains. Co-recorded alternatives within an archive × study × habitat panel do not imply local simultaneous effective pollination.

| Scope | panels | source rows | candidates with records | references with records | aquatic-positive C/R | mean recorded aquatic share C/R | same-panel other visitor among aquatic-positive C/R |
|---|---:|---:|---:|---:|---|---|---|
| All | 165 | 3837 | 25 | 76 | 10/33 | 0.045283 / 0.054434 | 10/10 ; 27/33 |
| Wet only | 7 | 327 | 9 | 17 | 5/10 | 0.113900 / 0.128537 | 5/5 ; 10/10 |
| All informative | 36 | 2045 | 21 | 51 | 9/28 | 0.056352 / 0.100150 | 9/9 ; 23/28 |
| Wet informative | 2 | 307 | 9 | 15 | 5/10 | 0.117172 / 0.147121 | 5/5 ; 10/10 |

Informative = at least one flagged visitor and both plant groups recorded within the same archive/study/habitat panel. Plants repeated across panels are aggregated before inference. Strict strata = region × island band × site band × recorded visitor richness (1, 2–4, 5–15, 16+); relaxed = region × richness. Mixed-label strata only; candidate-count-weighted mean within-stratum difference, 9999 label permutations with fixed LCG seed 7125. P-values nominal, exploratory, no multiplicity correction.

| Scope | Strict share difference | p | strict C/R | Relaxed share difference | p | relaxed C/R |
|---|---:|---:|---|---:|---:|---|
| All | -0.065156 | 0.0472 | 16/19 | -0.004868 | 0.8597 | 23/52 |
| Wet | -0.011534 | 0.8319 | 5/5 | -0.029430 | 0.6958 | 9/12 |
| All informative | -0.088367 | 0.0483 | 13/19 | -0.040376 | 0.3001 | 20/33 |
| Wet informative | -0.076818 | 0.4833 | 5/4 | -0.061841 | 0.4862 | 9/10 |

All directions are negative; nominal strict all-habitat signals are sensitive to stratification and are not confirmatory. No support for greater documented aquatic-hoverfly representation among wet candidates. Wet records are sparse and concentrated; this is not evidence of absence of aquatic-insect dependence.

Wet candidate positive taxa: Epilobium hirsutum (7/34 recorded visitor names flagged), Lythrum salicaria (2/17), Mentha aquatica (8/58), Rorippa amphibia (4/11), Crepis paludosa (1/5). All five have other visitors recorded in a panel that includes an aquatic-hoverfly record. Actual visitation rates, conspecific pollen contribution, selfing and replacement response are unknown. A functional dependence or light × aquatic dependency × redundancy model is not identified. Existing Page functional data and Lu illumination data have insufficient same-plant overlap; do not borrow guild-wide efficiency across plants. Existing macro result and main manuscript unchanged.

Replay from a repository with the frozen commit available:
```sh
node scripts/analysis/alan_wet_aquatic_visitor_support_20261124.mjs
```
