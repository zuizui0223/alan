# Official CLC-Change wetland vegetation sensitivity (source-only)

2026-11-20 | Original 17 CORINE 2000 initially mapped inland wetland/peatland botanical coordinates (2000 class 411/412).

Unlike comparing coarse status layers (minimum mapping unit 25 ha), **all three official 5-ha CORINE change epochs** are required: 2000–06, 2006–12, 2012–18. The official Copernicus FAQ recommends summing real transitions rather than subtracting 2000 and 2018 land-cover status inventories: https://land.copernicus.eu/en/faq/products/corine-land-cover .

Source endpoints (official EEA ArcGIS MapServer polygon layer 0):
- https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CHA2000_2006_WM/MapServer/0
- https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CHA2006_2012_WM/MapServer/0
- https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CHA2012_2018_WM/MapServer/0

For each of 17 pinned original inland wetland survey points, intersect actual original change polygons (code-before to code-after) with the **1 km ETRS89-LAEA circle** used as study-area approximation. Classify class 411 inland marsh and 412 peat bog gains and losses; other land-cover transitions are retained separately. Source version and exact polygon responses are SHA-digested. Neither the 25-ha status map nor the 5-ha change map detect all small wetland or vegetation changes. The source observation is a polygon-based spatial proxy, not causal exposure at a plant's exact coordinates.

**Keep independent from the original 187 matched-pair ALAN slopes.** The 187-pair source does not contain the complete original 374 survey-site coordinates. The prior original NDVI audit found 0 paired sites with complete vegetation histories. This 17-site mapping is an authentic supplementary habitat-change measurement, not an estimate on the original 187 pairs and not a substitute for a 374-site source recovery.

Reproduction: run the pinned GitHub Actions workflow and publish counts/areas only after ALL three official sources return valid polygons. The script throws an error and refuses an inferred all-zero result if the server schema, requested coverage, or geometries fail.


## Recovering the 374 original site coordinates

Added a new original-source recovery program: scripts/analysis/alan_recover_frozen187_pair_coordinates_20261120.py . It reconstructs EXACT 187 unchanged match identities from the pinned 1707-site source and joins only Git-SHA-pinned published coordinate tables: wet_annual/habitats, wet_shore/nearest_wells and seven_moth/nearest_wells. The program checks conflicting unit IDs; it prints aggregate genuine arm coverage, never approximates coordinates from the 100-km spatial block and never silently treats missing geography as (0,0).

An optional --original-geo-db path consumes the authoritative local original gbif_annual/annual.duckdb /geo table, if available, and certifies 187/187 fully recovered pairs. This source is not committed publicly. The original full-187 inference remains blocked until all 374 coordinate rows and matched-arm vegetation histories are genuinely recovered.
