# ALAN — Global same-island moth-record composition versus historically preceding artificial night light

2026-11-10. This is a source-verified OBSERVATIONAL detection/composition analysis; not insect density or nocturnal pollination.

## Source: actual 115 million original GBIF records processed, not sampling invented

- Frozen GBIF 2026-08 experimental DWCA source DOI: https://doi.org/10.15468/dl.8yrbe7. The four exact original Parquet families Noctuidae, Geometridae, Nymphalidae, Pieridae account for 115,233,537 raw rows of all years and basis of record. Only 2013–2025 georeferenced HUMAN_OBSERVATION with PRESENT status assigned to GSHHG 2.3.7 L1 islands >=20km2 entered the original sampled observational frame.
- Full original source extraction: https://github.com/zuizui0223/alan/actions/runs/42782146115 (SUCCESS), 22,713,063 geolocated source records, 1,217 genuine GSHHG independent islands, 228,361 distinct island × original 0.05-degree photo-site × observation-year cells; exact gzipped source SHA256 a993c0908600211e3e0ff66b88920ad0d47124e122f373bbce0e78ba9902ee6d.
- Within same photo-grid cell and year, 82,642 record-bearing cells have at least one of the two moth-family and two butterfly-family sources; 70,219 cell-years on 517 islands have both categories and >=5 total source records. These are photographic source records, NOT exactly 82,642 observed moths or actual moth–butterfly visits to flowers.
- GBIF source family counts in island polygons: Noctuidae 6,557,986; Geometridae 4,347,620; Nymphalidae 8,061,627; Pieridae 3,745,830 (sum 22,713,063). Other moth families are excluded rather than guessed. GSHHG 3,696 islands >=20km2 differs from GID 4,663 islands >=20km2; no silent identity crosswalk.

## New historical exposure from true original annual satellite pixels

Each real original island–cell–observation-year row is assigned an observation-year-minus-one corrected PCNL satellite pixel. The genuine source Li–Cao v3 annual 1992–2024 product DOI 10.5281/zenodo.17013414 is used for 2012–2024 images for observation years 2013–2025, respectively. All 13 source rasters have frozen original MD5 values; each worker refuses source changes.

The original GBIF records were assigned to islands by each true coordinate, but the published large-scale source table stores only the 0.05-degree grid ID. Therefore sampling the PCNL value at the grid CENTER is an approximation and the center MUST ALSO be inside the SAME original GSHHG land polygon. An original cell whose center falls into ocean is NA; it is not a dark source site, not an absent moth, and not a source error. This restriction will disadvantage narrow/coastal islands; report original and retained site/island coverage and area sensitivity. No 2024 VIIRS photo light is backprojected onto 2013 original insect records.

Every matching site retains original true island ID, source observation year, source 0.05-degree cell indices, original moth/butterfly-family record counts, satellite history-year and true numerically sampled pixel value. The collector demands all 13 source raster years across four year shards, exact input SHA, original 228,361 row keys and original 22,713,063 observational record sum. Missing pixel remains NaN, not zero.

## Exploratory true global observational contrast

Primary estimand: (Noctuidae + Geometridae) GBIF source record fraction among these two and (Nymphalidae + Pieridae) in a REAL same-island, same-observation-year source photo cell with BOTH categories documented and >=5 original source observation records. This is an ascertainment/composition estimand, NOT a pollination preference, true local presence/absence, insect biological population trend or absolute abundance.

Observational model: source moth-family PHOTO-RECORD percentage ~ log1p(previous-year real PCNL satellite light) + log1p(total four-family genuine GBIF documentation opportunities) + source site-center latitude; within-island × observed-year intercepts, equal total weight per actual island and nominal island-clustered covariance. Require >=2 distinct real 0.05-degree geographic cells in that island/year and >=0.02 source loglight span. No fitting just to create a significant p-value.

Sensitivities: each photo-cell year equal, each original island-year equal, without total source-effort adjustment, include source-zero category (recorded no moth is still not true biologic absence), >=20 documented events, island-area × PCNL interaction, largest-observed-island exclusion. Neg control: light slope against total four-family GBIF source-record count itself, which can reveal brightness-associated photography effort.

Do not treat a positive/negative record fraction coefficient as caused by moth attraction, altered actual moth density, seed set, pollen deposition or floral evolution. The source is uncontrolled for dataset-specific platform, observer, capture availability, horticulture, month/season and urban land-use, and may not be representative of plant-associated moths. A mainland comparison must use independent true mainland source rows (NOT infer it from the island-only source or a country total). A simultaneous real moth-monitoring outcome needs separate standardized censuses (e.g. Korean island field dataset) and still controls phototaxis.

## Publication gate

Only report an ecological signal if after exact source-land centers, true year-prior raster joins, source observation effort control, island-year comparisons, and independent/small-island sensitivity the direction remains stable. Even this does not prove causality; confirm actual insect activity and effective floral pollination independently. If source support disappears or coefficients depend on source visibility, report the observational failure mode as the main result rather than upgrading GloBI photo source labels.

Global source year-PCNL workflow: https://github.com/zuizui0223/alan/actions/runs/42786502922. Source reference is complete only if all four worker shards and collector report success.