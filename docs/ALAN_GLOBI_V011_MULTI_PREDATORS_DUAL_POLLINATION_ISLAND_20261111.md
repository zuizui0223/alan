# ALAN / full GloBI v0.11 — predator communities and multiple ecological roles

**2026-11-11: actual complete 24,453,954-row GloBI source executed.**
Research branch, not a confirmed artificial-nightlight impact on pollination or plant loss.

## Main empirical finding: a predator can also be a flower visitor

Existing independently completed signed GitHub Actions run:
https://github.com/zuizui0223/alan/actions/runs/42944244963

Release: https://zenodo.org/records/22691479
Source SHA256: 670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3

| Guild | Animal species names with family-confirmed moth predation/kill relations | Flower visit or database pollination species | SAME named species with broad moth-or-unresolved-Lepidoptera predation and flower use |
|---|---:|---:|---:|
| Bats | 94 | 55 | 17 |
| Birds | 349 | 505 | 157 |
| Geckos | 7 | 4 | 1 |
| Centipedes | 3 | 2 | 1 |
| Spiders | 155 | 90 | 40 |
| Mantises | 17 | 10 | 6 |
| Predatory wasps | 58 | 644 | 88 |
| Frogs | 3 | 2 | 0 |

This **broad source screen includes** eats, preysOn, kills and hasPrey, and the dual-role count also includes prey classified to Lepidoptera without a resolved family. Thus these counts are not proof of direct predation on adult flower-visiting moths.

Source details expose why the functional qualification matters:
- Hemidactylus frenatus has moth predation records plus a GloBI database-pollinates label, but the latter requires original-study verification: nectar intake is not proven effective pollen deposition.
- Phelsuma ornata and P. cepediana have flower visits and reported pollination; not all geckos should be assigned the same foraging role. Independent field/experimental studies demonstrate gecko nectar use.
- Thereuopoda clunifera is the recorded dual-role centipede, but the predation link is family-unresolved Lepidoptera and its flower link is a VISIT. GloBI reports ZERO centipede species with a database pollinates relationship. It is not acceptable to claim centipedes are effective plant pollinators.
- Even spider flower visits can be flowers used to hunt prey rather than pollen transfer.

The new strict reproducible GloBI version-locked scanner in this branch distinguishes eats/preysOn from kills, butterfly vs taxonomically identified moth, reported effective pollination vs visits, plant dispersal, and animal-to-animal predation among predator guilds, retaining source and life-stage limitations:
scripts/analysis/alan_globi_predator_pollinator_multilayer_v011.py
tests/test_alan_globi_predator_pollinator_multilayer_v011.py
.github/workflows/alan-globi-v011-predator-pollinator-roles.yml

## Independently executed: exact same-island gecko + full GloBI plants + photos + VIIRS

- Official source WCVP crosswalk of the 1,651 positive GloBI v0.11 moth-visited plant names matched **1,356** accepted name entries.
- GBIF original 2019/2024 HUMAN observations of Hemidactylus frenatus and H. mabouia, <=1 km positional accuracy, assigned to the independently pinned original GSHHG 2.3.7 landmass ID.
- Plants from one-photo-per-observation original historical_main source, >=3 authentic photo opportunities, >=2 accepted photographed taxa, real 2019/2024 same-year VIIRS at least three valid sampled months.
- Geographic overlap requires the SAME year, the SAME physical island landmass ID and the SAME 1-degree lat/lon grid. Ambiguous 0.05-degree GSHHG landmass cells excluded.

| Source year | Land class | Qualified photo site-years | Sites with same island/geographic positive GBIF gecko records | Of those with photographed GloBI moth-positive plant species | GSHHG landmasses with both |
|---|---|---:|---:|---:|---:|
| 2019 | Island | 505 | 169 | **52** | 26 |
| 2019 | Mainland | 2,714 | 271 | **89** | 5 |
| 2024 | Island | 1,349 | 348 | **107** | 55 |
| 2024 | Mainland | 8,793 | 1,454 | **432** | 5 |

The derived botanical observation-based co-location is stronger evidence than a global moth-vistor taxon list alone, but is NOT a local gecko predation event, and the coarse within-island 1° grid can still separate the animals by tens of kilometers. A negative GBIF record must never be inferred from a missing record. No causal mediation, gecko abundance or annual VNL ALAN model is estimated.

Machine-readable public aggregate (no private plant/photo site rows):
results/alan_globi_full_v011_gecko_same_island_photo_source_20261111.json

## Five linked but non-equivalent functional layers

1. Predation on moths: actual documented eats/preysOn, separating known adult moths from larvae/unspecified prey.
2. Flower attendance by predators: visitsFlowersOf, not automatically pollination.
3. Pollination claims: GloBI pollinates relationship, followed by independent efficacy/seed-set evidence before assigning functional service.
4. Seed dispersal by predators such as fruit/nectar bats or geckos: a distinct relationship.
5. Predation ON predators: centipedes sometimes take geckos or bats, allowing an intraguild cascade.

The ecological hypothesis is conditional, not a discovery: ALAN may redistribute nocturnal moths and their predators, and both capture risk and alternative pollination functions may differ on islands vs mainland. Positive-only original GBIF predator records are useful to locate where to inspect this mechanism; the lack of records is not confirmed predator absence.

## External sources for interpretation

- ALAN and nocturnal moth pollination review: https://doi.org/10.1111/een.12174
- Centipede predation on bats/reptiles review: https://doi.org/10.3389/fevo.2025.1634037
- Independent centipede feeding source with documented insect predation and some plant material consumption, but no proof of flower pollination: https://pmc.ncbi.nlm.nih.gov/articles/PMC12782032/
- Phelsuma ornata experimental flower-resource choice: https://pmc.ncbi.nlm.nih.gov/articles/PMC1618895/
- Hemidactylus frenatus nocturnal hunting and incidental nectar consumption: https://www.iucngisd.org/gisd/speciesname/Hemidactylus%2Bfrenatus

Avoid significance hunting and mixing global GloBI relationships with unobserved site-specific predation, native/invasive status, or effective pollen deposition. Current research findings are source coverage and functional network opportunity, not an island ALAN causal coefficient. Preserve Draft status.


## 2026-11-11 completed strict GloBI source: stronger than initial broad counts

**Completed SHA-locked processing, independently verified**: The stricter source runner scanned exactly all **24,453,954** original v0.11 interactions; source digest and all four output members match. Full-source workflow run [42952223131](https://github.com/zuizui0223/alan/actions/runs/42952223131) ended *FAILURE only after scan and artifact upload*, because its final test referenced an incorrectly named JSON caveat key. The artifact **13826322051** was produced intact; a separate no-redownload verification run [42954755175](https://github.com/zuizui0223/alan/actions/runs/42954755175) completed **SUCCESS**, SHA-verifying each extracted file and confirming source row totals and caveats. The original workflow is now corrected and manual-only, avoiding another costly 2.5GB download.

The stricter scan excludes kill-only claims from its moth-predation arm, separates butterfly prey, records moth prey life stage, and distinguishes flower visits from database-pollinates claims.

| Animal group | Source-positive species with moth predation | Source-positive species with flower visit or DB pollinates | SAME source species with both functions | Both after iNaturalist+trophiCH exclusion |
|---|---:|---:|---:|---:|
| Bats | 94 | 54 | **6** | **6** |
| Geckos | 7 | 4 | 1 | 0 |
| Centipedes | 3 | 2 | **0** | 0 |
| Spiders | 155 | 75 | 30 | 1 |
| Birds | 349 | 435 | 84 | 3 |
| Mantises | 17 | 10 | 6 | 0 |

Among **41,270** source moth predation relationship rows, documented prey stage is **adult 4,420**, **nonadult 547**, **unknown 36,303**. Therefore the large majority of source claims do **not** document killing flower-visiting adult moths. The strict 6 bats are species names with separate family-confirmed moth prey and flower interactions (in this dataset they also have a DB-pollinates relationship), NOT confirmed adult moth capture and effective local flower fertilization.

Positive plant-name associations by animal group: bats **135** plant names, geckos **10**, centipedes **1**, spiders **167**, birds **1,170**. For centipedes, **zero DB-pollinates** relationships; a flower visit is not proven service.

Public summarized strict source data and exact output-member hashes:
[results/alan_GloBI_v011_strict_predator_pollinator_real_all_24m_20261111.json](https://github.com/zuizui0223/alan/blob/research/moth-gecko-flora-feasibility-20261111/results/alan_GloBI_v011_strict_predator_pollinator_real_all_24m_20261111.json).

## First actual non-gecko GBIF geographic observation cohort

Successful [real GBIF source run 42954162844](https://github.com/zuizui0223/alan/actions/runs/42954162844), original aggregate-only artifact **13826800146**, read 2019 and 2024 exact direct positive observations with positional uncertainty <=25km:

| Target animal taxon | 2019 recorded directly | 2024 recorded directly |
|---|---:|---:|
| Glossophaga soricina (bat) | 51 | 38 |
| Antrozous pallidus (bat) | 56 | 82 |
| Mystacina tuberculata (bat) | 3 | 3 |
| Thereuopoda clunifera (centipede) | 6 | 53 |
| Scolopendra subspinipes (centipede) | 62 | 173 |

Vespa velutina, Argiope bruennichi and Mantis religiosa cannot be fully counted under the bounded 12k source-event cap; all remain **UNRESOLVED**, NOT "zero" occurrence. Data are public source observations but this study deliberately does not export rare bat exact coordinates. Their growth is NOT an organism abundance trend; GBIF monitoring effort is uncontrolled. Source code at scripts/acquisition/alan_multiguild_gbif_positive_pilot.py and tests/test_alan_multiguild_gbif_positive_pilot.py. Raw GBIF positive records are NOT observed predation or visitation events.

**Completed source stage, not biological mechanism**: no valid predator detection–nondetection denominator, same-plant predator+visitor interaction or long-term annual light exposure has yet been observed across all predator groups. The next model requires same-landmass/polygon geography plus reliable observation effort and preferably verified adult moth capture vs flowers' pollen deposition. No predator-group ALAN regression coefficient or p value is claimed here.
