# Site context and moderator readiness

Updated2026-11-20. `primary_coordinates.csv` preserves primary-source site locations, precision/extent, year, experimental setting and unresolved limitations. `catalogue_coordinate_coverage.csv` joins these to current programme labels without multiplying effect rows by sites.

Seven of the ten catalogue programme labels now have at least one numeric site coordinate. This is not seven approved independent studies or seven remotely measured wetness values. Giavi contributes12 meadow locations within one programme; the other six located programmes have representative site coordinates. Guo has four named sites but no verified numeric coordinates here; Hey has a named farm; Sanders2016 is still absent from this coordinate table. Missing coordinates do not mean dry habitat or ineligibility.

## Newly verified sources

- Bennie2018, [primary paper](https://doi.org/10.1111/1365-2664.12927), Methods2.1; cached `grassland2018/primary.txt`:50°2′N5°12′W. The35×35m field enclosure has shallow poorly drained gley soil. Coordinate precision is one arcminute, much coarser than the experiment. Do not classify this grassland as dry or assign plot-nearest-water distances from that point.
- Parkinson2023, [primary Methods2a–b](https://pmc.ncbi.nlm.nih.gov/articles/PMC10613537/): Waterfowl Pond42°55′51.2″N83°11′23.9″W; overnight3–4June2022. One coordinate represents adjacent littoral and riparian treatments, not each trap. Their contrast is not a general wet-versus-dry comparison. The cached `parkinson2023/article.xml` is a JSON500 error response, not the article; this update used the readable primary PMC page.

Species EIV MOIST remains separate from site wetness. Existing Giavi GSW values describe1984–2024 background water occurrence, not2016 conditions; they cannot serve as same-year hydrology. No new remote sensing values, climatic moderators, or moisture meta-regression are produced in this update. Finite species-moisture support remains10 effects/5species/3programmes, as documented in `plant_moisture_manifest.json`.

The coordinate coverage table is an audit of current inputs, not a new effect-size catalogue. Each listed source location and unresolved field must be reviewed before remote covariates are attached; treatment locations or sufficiently defined experimental extents are required for fine-scale exposure interpretation.

## Annual climate context acquired (2026-11-20)

`annual_climate.csv` now contains25 site-years at18 named locations in7 catalogue programmes;15 distinct returned grid-year combinations. `daily_climate.csv.gz` retains the daily inputs. This supersedes the earlier no-climate-acquisition statement above, but not the wetness limitations.

[Open-Meteo historical API](https://open-meteo.com/en/docs/historical-weather-api) was explicitly queried for ERA5 (0.25-degree), GMT daily mean temperature and total precipitation. Nearest-grid selection and disabled elevation downscaling preserve grid context. Annual means/sums require complete365/366-day records; no missing days were filled. Returned grid coordinates, distance from requested points, request URLs, retrieval times and raw hashes are retained in `climate_requests.csv` and `climate_audit.json`.

Calendar years follow the primary coordinate table, including every year in2012–2016 or2012–2015 ranges. They are not exact experimental exposure windows, long-term climate normals or pre-treatment baselines. Named locations and shared grid-years are not additional independent experiments. Original coordinate precision/datum uncertainties persist. Artificial mesocosm watering, substrate, lake water levels and local inundation are not inferred from grid rainfall.

Reproduce with `python database/alan_experiments/site_context/acquire_annual_climate.py`; cached raw responses are reused. No effect-level averaging across years/sites, wet/dry classification or meta-regression is performed. Hey/Guo and absent Sanders coordinates remain unfilled. Source attribution: ERA5 reanalysis distributed through Open-Meteo; API documentation linked above. This adds climate context, not a completed environmental-moderator analysis.

Independent agent verification matched18 raw hashes, all9,142 daily records and25 annual aggregates. Eight site-years contain365days and17 contain366days. Maximum recomputation differences were1.78e-15 degrees C and1.36e-12mm. Source coordinates, requested years, returned grids, units and grid distances also matched. This is computational validation, not independent meteorological observation.

## Within-pair spatial support

`climate_pair_support.csv` audits all six Giavi2016 meadow pairs against the returned ERA5 grids. Five pairs share one grid cell within each pair: their identical climate values do not establish identical local conditions. Pair3 spans two cells. Its requested points are1.810km apart, while returned cell centres are27.799km apart. The annual lit-minus-dark differences are-3.898 degrees C and+328.5mm. These are grid-context contrasts, not measured plot differences or treatment effects. Real topographic climatic variation may exist, but these data cannot resolve the within-pair contrast.

No within-pair climate adjustment is adopted from these values. No point is deleted, snapped to its paired point, or assigned a corrected value. The annual context remains useful at its stated coarse support; it cannot establish climate matching of the illumination/control plots. Reproduce with `python database/alan_experiments/site_context/audit_climate_pair_support.py`; the audit hashes its annual source. A later finer-scale source requires explicit spatial and temporal validation rather than selection for agreement with the light response.

Independent agent recalculation matched all six pair assignments, differences, distances and the annual-input hash. This verifies the spatial-support audit, not local meteorological accuracy.
