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


## Executed source-authentic results (GitHub Actions #44668673040, SUCCESS)

The three officially published EEA polygon sources all returned valid GeoJSON, and five unit tests of exact transition and circular-area calculations and four unit tests of source-site identities completed successfully. Full original 51 site-epoch records, actual transition codes and per-request response SHA256 digests are preserved at [CORINE result](../results/all_species_forest/wet_falsification/corine_official_cha_17wetland_20261120.json); [original-coordinate support](../results/all_species_forest/wet_falsification/frozen187_original_coordinate_recovery_20261120.json).

- 17 originally mapped (2000) inland marsh/peatland botanical points × three official change periods = **51 EEA change-epoch queries**.
- **Zero** sites with mapped inland-wetland 411/412 loss within their 1-km circles; **one** site with mapped wetland-class gain, **60.52 hectares** (ID 4297, 2012–2018 mixed forest 313 → peat bog 412).
- Four sites have any mapped land-cover change: 4297 above; 2191 had coniferous forest 312→transitional woodland/shrub 324 (5.26 ha, 2012–2018); 7816 had transitional woodland/shrub 324→moors/heathland 322 (3.79 ha, 2012–2018); 4685 had complex cultivation 242→industrial/commercial 121 in 2000–2006 (6.87 ha) and 2006–2012 (9.31 ha).
- Land-cover transitions at nonwetland locations of the 1-km circle are NOT wetland 411/412 gain/loss. These are spatially clipped areas, not a mapped area measured on the individual plants.

### Source recovery for same 187 matched pairs

From three existing source-hash-verified site-coordinate tables there are **109 unique recoverable survey coordinates overall**; only **3 brightening arm** and **1 comparison arm** coordinates belong to the original 374 distinct matched sites. There are **zero complete geographic matched pairs**. The authoritative full geo table was NOT present in the repository; original local gbif_annual/annual.duckdb is required. Avoid manufacturing a 187-pair vegetation-adjusted coefficient from these 17 independent wetland sites.

**Interpretation**: No CORINE-mapped ≥5-ha inland-wetland loss was observed near these seventeen locations. This weakens only the *coarse mapped loss* explanation for these locations; it cannot exclude changes of wet grassland, drainage, groundwater, 5-ha-subthreshold vegetation, land-management history or correlated survey effort. No light-specific causal coefficient is identified.
