# ALAN x fd x island — moth-like phenotype, geckos and artificial nightlight

**Status: runnable source-support/validation scripts, not an ecological effect.**
Independent ALAN research branch: read the frozen public fd colour/form reference and authorized island identity-level sources, without changing island Chapter 1 or earlier ALAN white-ROI results.

## Script 1: morphology and GloBI

Run:

    python scripts/analysis/alan_moth_syndrome_validate.py --profiles /secure/profile_coordinates.csv --traits /secure/reviewed_taxon_traits.csv.gz --interactions /secure/GloBI_bounded_interaction_source_cache.csv.gz --output /secure/moth_validation

Inputs:
- profiles: public fd upstream floral_space_20261106/all/profile_coordinates.csv, columns profile_id,flower_primary_color,floral_form (pipe-separated state sets).
- traits: authorized local species/direct evidence, accepted_species,flower_primary_color,floral_form,tube_depth_class. This taxon-level table cannot be inverted from public fd island×profile counts. Resolve names and verify evidence before use. Do NOT commit raw data.
- interactions: positive source-frozen GloBI plant_name,guild,relation records (the previous ALAN source cache can work). Existing GloBI queries were bounded; a record is positive flower visitation, not effective nocturnal pollination or proof of moth dependence.

The strict candidate is WHITE plus TUBULAR/SALVERFORM/FUNNEL_TRUMPET plus DEEP tube. A multistate or missing trait stays UNKNOWN, not a negative. No Cartesian product of recorded colour/form sets is generated; green_brown_inconspicuous is not automatically pale. The output is an aggregate validation table and a 496-profile colour/form PREFILTER, which does not contain deep-tube information.

Butterfly-only recorded does not imply moth absence; taxa with both GloBI guilds remain separate. Descriptive odds ratios are NOT biological specificity.

## Script 2: source overlap

Run:

    python scripts/analysis/alan_island_gecko_light_support.py --profiles /secure/profile_coordinates.csv --traits /secure/reviewed_taxon_traits.csv.gz --floras /secure/reviewed_island_flora.csv.gz --geckos /secure/reviewed_island_gecko.csv --light /secure/island_annual_viirs.csv --islands /secure/island_metadata.csv --output /secure/gecko_support

Additional columns:
- floras: island_id,accepted_species,establishment_status,review_status[,establishment_year]. Island-specific reviewed naturalized means naturalized. Cultivated/transient/unknown are excluded, and mere introduction is not proof of naturalization. First GBIF observation cannot substitute for documented establishment year.
- geckos: island_id,gecko_status,review_status,geographic_scope,source. Accept only reviewed island-scope present/confirmed_absent/unknown. Missing observations and country GRIIS checklists cannot demonstrate absence or presence on individual islands. Native AND introduced geckos may matter.
- light: island_id,year,radiance,product_version,n_valid_months,coverage_fraction. Use actual annual VIIRS radiance measured over each island polygon, default 2019–2025, minimum 6 months, 70% area/pixel coverage and two qualified years. Do not substitute quarterly source-photo-cell light.
- islands: island_id,archipelago,area_km2 from pinned GSHHG geography. Join by unique polygon IDs.

Output: support_aggregate.json and PRIVATE_island_derived_counts_DO_NOT_COMMIT.csv. Keep the latter LOCAL, not in a public repository: it merges identity-level derived restricted trait evidence. Publication of even aggregates requires license and small-cell checks.

The 2×2 gecko × median-nightlight table is a SOURCE SUPPORT diagnostic ONLY. The subsequent ecological model, only if source size and independently represented archipelagos are adequate, should use continuous log(1+radiance), moth-like trait × gecko × nightlight, spatial and taxonomic clustering and covariates for development/ports, introduction pressure, native floristic pool, area, isolation and trait recording effort. True independent daylight flower-visitor evidence is a better negative control than non-moth-like shape.

**Temporal caveat:** Present naturalization may predate all modern lighting. A present-day checklist versus contemporary satellite light yields cross-sectional association, not colonization failure due to light. Separately count documented establishment since 2013 (not equivalent to illumination onset). DMSP-OLS begins in 1992, but is not directly radiometrically comparable to VIIRS; explicit calibration is required for a combined trend.

No ecological coefficient, significance or complete-case island count has yet been computed: reviewed island-level gecko absence, truly naturalized flora and aligned annual nightlights must be supplied first. The original ALAN photo-colour ROI remains an entirely different response.

Tests:

    python -m unittest discover -s tests -p 'test_alan_moth_gecko_flora_support.py' -v
