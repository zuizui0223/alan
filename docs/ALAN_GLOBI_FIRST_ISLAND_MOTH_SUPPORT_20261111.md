# GloBI-FIRST island moth–plant × nightlight: real-source audit (2026-11-11)

**Status: real archived datasets joined and audited, not a fitted ALAN/gecko effect.** This is an independent research branch of ALAN. No changes were made to fd or island Chapter 1.

## Immutable sources

- mois Actions run 41976552977, artifact 13391803415: 61,967 actual moth/butterfly interaction claims from 7 moth families and 5 butterfly families, source SHA a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994.
- island Actions run 34232450884, artifact 10058653212: 1,039,757 historical island × plant-name rows, source SHA 6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a.
- mois Actions run 42171173575, artifact 13481850056: 889 photo islands, 3,575 island-years, 61,734 original plant photos, source SHA acaf9287a29ef01e2fc79fd4b3450b4c09bb5eee16737e55b4326a80359cbdf8.

The new source reader fails closed on SHA, schema, duplicate island–taxon rows, unexpected origin labels or incorrect VIIRS-product claims. Only source-supported aggregate counts are published, never the private original species×island ledger.

## Executed actual-source findings

| Source evidence | GloBI moth-associated plant names | Exact matches in island floras | Islands with >=1 matched moth-positive plant | Islands with >=5 |
|---|---:|---:|---:|---:|
| All stage-2 cached source records | 2,362 | 1,880 | 3,404 | 2,417 |
| Exclude iNaturalist interaction source | 1,057 | 821 | 2,890 | 1,881 |
| Exclude both iNaturalist and trophiCH interaction sources | 786 | 634 | 2,850 | 1,848 |
| Exclude iNaturalist and restrict to visitsFlowersOf records | 540 | 469 | 2,679 | 1,728 |

For the all-source match, the 1,880 species names occur on 175,485 historical island×name rows. Origin status: 13,433 native, 4,745 introduced, 157,307 unresolved. Introduced is NOT naturalized; 688 moth-visit-associated plant names have introduced records on 121 islands.

Among the total GloBI cache, 1,843 plant names have both recorded moth and butterfly visits; 519 have moth but no *recorded* butterfly visits. Recorded-only is NOT biological exclusivity. The absolute number of interaction rows (13,262 moth / 48,705 butterfly) is not a count of independent experiments, a count of species, or a measure of interaction strength.

Of the 61,734 original plant photographs, 14,041 show plants whose species have a positive moth visit record SOMEWHERE in GloBI (NOT a photo depicting a moth visit). These occur on 590 islands; 213 islands have >=5 such plant photos, 72 have >=20, and 344 have >=2 photo years. The existing same-year light data are four-month quarterly proxies, not annual VIIRS: 19 photo islands / 25 island-years, with moth-tagged plant photos on 24 of those island-years.

All previously published fd flower-shape prefilters are now optional exploratory modifiers, NOT entry criteria or proof of moth pollination.

## What the findings do NOT demonstrate

A plant with a documented global moth visitor is not proven to use moths on each of its islands. Unrecorded interactions are not biological absences. A GloBI interaction derived from iNaturalist may not be independent of another iNaturalist plant photo, which is why source-exclusion sensitivity matters. GloBI pollinates is also NOT automatically direct nocturnal stigma pollen deposition or seed set. Historical island introduced presence is not naturalization and a cross-sectional species list cannot document recent ALAN-caused extinction.

No artificial-light × gecko × moth-pollination coefficient or p value was fitted.

## Full GloBI stable release: next large acquisition

GloBI explicitly advises research projects to use stable versioned products rather than unversioned live APIs: https://www.globalbioticinteractions.org/data .

Verified September 10 2026 interpreted release v0.11: https://zenodo.org/records/22691479 ; DOI 10.5281/zenodo.22691479 . The complete interactions.csv.gz is about **2.5 GB compressed**, with MD5 df4a2adaef38b848d00556df5019ed3e and SHA256 670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3.

**It was NOT downloaded or processed for the present source audit.** The count of 3,404 islands refers only to the older bounded 61,967-row live API cache (seven moth families and five butterfly families). For a final source run, stream and verify the stable version, evaluate full lepidopteran taxonomic coverage, preserve study provenance and interaction event geolocation/date when actually supplied, and crosswalk WCVP accepted names rather than assuming all unrecognized names are biological absences.

## Scientific scale and next empirical estimand

- **Background:** island plant lists × global GloBI-positive interaction evidence describes a long-term species-pool pattern. It cannot be used as an annual light response.
- **Time-matched response:** original accepted plant species photographed or recorded in the same island/site/year as genuinely annual VIIRS, with moth-positive AND butterfly-positive flags, overlap retained, observation effort, region and community composition controlled. Report as observed photographic plant-assemblage composition (not local actual moth pollination).
- **Mechanism:** genuine georeferenced, dated interaction observations, nocturnal pollen transfer or effective pollination, and independent gecko surveys with confirmed absence. This is not yet available as a complete join.
- **Separate naturalization:** exact GloNAF naturalized floras and geographically audited island regions. Do not promote historical introduced flags to naturalized.

This branch only improves source support and enables an independent GloBI-first macro research design. No positive light-mediated ecological mechanism has yet been inferred.

## Reproduction

Run the new source script with the three frozen artifact ZIPs:

    python scripts/analysis/alan_globi_primary_island_support.py --globi-zip SOURCE_GLOBI.zip --island-zip SOURCE_ISLAND.zip --photo-zip SOURCE_PHOTO.zip --output output.json

Five local synthetic tests of the equivalent source-join logic passed, and the actual source-joined aggregate run completed. See results/alan_globi_primary_island_source_verified_20261111.json for published counts.
