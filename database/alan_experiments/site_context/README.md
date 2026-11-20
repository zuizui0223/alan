# Site context and moderator readiness

Updated2026-11-20. `primary_coordinates.csv` preserves primary-source site locations, precision/extent, year, experimental setting and unresolved limitations. `catalogue_coordinate_coverage.csv` joins these to current programme labels without multiplying effect rows by sites.

Seven of the ten catalogue programme labels now have at least one numeric site coordinate. This is not seven approved independent studies or seven remotely measured wetness values. Giavi contributes12 meadow locations within one programme; the other six located programmes have representative site coordinates. Guo has four named sites but no verified numeric coordinates here; Hey has a named farm; Sanders2016 is still absent from this coordinate table. Missing coordinates do not mean dry habitat or ineligibility.

## Newly verified sources

- Bennie2018, [primary paper](https://doi.org/10.1111/1365-2664.12927), Methods2.1; cached `grassland2018/primary.txt`:50°2′N5°12′W. The35×35m field enclosure has shallow poorly drained gley soil. Coordinate precision is one arcminute, much coarser than the experiment. Do not classify this grassland as dry or assign plot-nearest-water distances from that point.
- Parkinson2023, [primary Methods2a–b](https://pmc.ncbi.nlm.nih.gov/articles/PMC10613537/): Waterfowl Pond42°55′51.2″N83°11′23.9″W; overnight3–4June2022. One coordinate represents adjacent littoral and riparian treatments, not each trap. Their contrast is not a general wet-versus-dry comparison. The cached `parkinson2023/article.xml` is a JSON500 error response, not the article; this update used the readable primary PMC page.

Species EIV MOIST remains separate from site wetness. Existing Giavi GSW values describe1984–2024 background water occurrence, not2016 conditions; they cannot serve as same-year hydrology. No new remote sensing values, climatic moderators, or moisture meta-regression are produced in this update. Finite species-moisture support remains10 effects/5species/3programmes, as documented in `plant_moisture_manifest.json`.

The coordinate coverage table is an audit of current inputs, not a new effect-size catalogue. Each listed source location and unresolved field must be reviewed before remote covariates are attached; treatment locations or sufficiently defined experimental extents are required for fine-scale exposure interpretation.
