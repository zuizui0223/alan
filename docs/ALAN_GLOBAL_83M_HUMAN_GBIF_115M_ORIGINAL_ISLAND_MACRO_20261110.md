# ALAN | Worldwide moth and butterfly observation macro-analysis: real-source expansion

Updated 2026-11-10. This is a real-data source and coverage program, not yet evidence of light-caused pollinator population change.

## Data and verification

GBIF public human-observation indexed counts, georeferenced, 2013–2025, source-verified in Actions #42780755320:
- 83,309,815 Lepidoptera occurrence records
- 41,838,694 butterfly records in six verified families
- 41,471,121 remaining Lepidoptera observations, not all verified nocturnal moths
- 36,436,081 observations in the independently resolved subset of selected moth families (Saturniidae ambiguous in live GBIF taxonomy and excluded rather than guessed)
- 138,688,351 georeferenced human Insecta occurrences (not an effort-controlled denominator).

All these are publication-indexed OBSERVATION records, not captured population numbers, effective floral pollen transfer or plant seed set. Country-level counts include mainland and cannot stand in for island support.

The frozen 2026-08-01 GBIF DWCA preview has DOI https://doi.org/10.15468/dl.8yrbe7, partitioned by taxon on public S3. Parquet metadata from the actual public files verified in Actions #42780706666: 170,319,717 Lepidoptera all-era/all-basis occurrence rows, 6 taxon partitions, total 39.330 GiB compressed. This is not the same denominator as the live 2013–2025 human-only GBIF API count.

Four genuine original taxonomic source partitions selected to validate observed guild-record contrast:
- Noctuidae 33,646,170 original rows (moth)
- Geometridae 24,782,143 original rows (moth)
- Nymphalidae 40,565,469 original rows (butterfly)
- Pieridae 16,239,755 original rows (butterfly)
These are 115,233,537 raw all-years/all-basis rows before original-date, record-type and location restrictions.

The first two original-row-group real-geography pilot (Actions #42781692909) read exactly 200,704 source rows per tested family, with actual original point-in-polygon GSHHG membership:
- Noctuidae: 60,037 authentic georeferenced human observations across 58 GSHHG islands, at 11,247 unique original 0.05-degree source site-year cells
- Nymphalidae: 35,841 authentic georeferenced human observations across 84 GSHHG islands, at 17,403 original 0.05-degree source site-year cells
- GSHHG 2.3.7 L1 landmasses with verified polygon area >=20 km2 number 3,696; this is not the independently sourced 4,663-island GID universe.

These source-row-group examples are intentionally sorted-data pilots; they are NOT representative geographic frequencies, independent datasets, or a final global estimate.

## Actual complete-source method

Stream all original Parquet row groups across the four disjoint taxa from a byte/version-pinned GBIF source, only actual 2013–2025 HUMAN_OBSERVATION records, PRESENT and georeferenced; identify original GSHHG polygon membership at each real observation point before aggregating to 0.05-degree source cell x observation year x island x original moth or butterfly family.

Store source counts and distinct taxonomic keys per exact island-cell-year without republishing raw contributor positions. Reject overlapping ambiguous island geometry and retain unknown taxonomic classifications. Count co-supported cells where two true family guilds have observed source records and at least five source events, but NEVER treat unrecorded species as biologically absent.

Then match actual observation year minus one to independently measured annual satellite nighttime light (2024 VIIRS cannot be used for observations in earlier calendar years), with urban land use, region, dataset, visitor family and observation effort controls. Compare islands and mainland at equal within-region/date support and cluster by island/region. Full original source processing must pass complete Parquet-row, source-key and group-collision checks before ecological modeling.

Critically the global GBIF photograph dataset is a nonstandardized citizen science detection dataset. A positive association might reflect light attracting moths or attracting OBSERVERS, and recording/horticulture shifts; it is NOT demonstrated population growth or effective nocturnal pollination. Independent actual light-trap studies (e.g. 30 Korean islands/30,442 captured moths) are separate validation and do not magically correct this recording bias.

## Source-snapshot caveat

GBIF issue #609 reports about 4.9 percent missing files in regular public AWS GBIF monthly snapshots since July 2025, including October 2026. For this experiment we use the separately derived 2026-08 DWCA preview, not the potentially incomplete routine 2026-10 AWS snapshot. The preview is experimental and requires independent full-family row-count and release provenance checks.

Original true global source workflow: https://github.com/zuizui0223/alan/actions/runs/42782146115. Full support values are only valid if all four family workers and the final collector complete source checks.