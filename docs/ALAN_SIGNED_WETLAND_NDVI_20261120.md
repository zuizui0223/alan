# Wetland vegetation changes not resolved by coarse CLC class maps

2026-11-20 | supplementary original-source analysis (NOT 187-pair adjustment).

The fixed 17 originally CORINE-2000 411/412 inland wetland/peatland
survey sites have **zero mapped wetland-class losses** in the official
2000–2006, 2006–2012 and 2012–2018 minimum-5-ha change products.
But lack of land-cover-class loss is not evidence that green biomass,
flooding, hydrology, or mowing has stayed constant.

Here independently recompute 2000–2023 Terra MOD13Q1 250-m, 16-day
growing-season QA0 NDVI. Require monthly pairs 4–5, 6–7, 8–9 and
>=12 full-season years across >=18 years, >=2 years from 2000–04 and
>=2 from 2019–23, matching existing original-source audit. Compute
signed OLS and robust Theil–Sen slope per decade plus actual early
versus late 5-year-window median differences.

Join the original unchanged 17 botanical site-level candidate-minus-
reference recorded-rate trends and the official CLC change-polygon
areas, both by genuine original unit ID. Give site counts and spatial
block bootstrap intervals; only 8 100km spatial blocks, postselected
candidate species, and one site of long-term satellite brightening
preclude causal ALAN interpretation.

Important scale separation: the MODIS vegetation signal is a 250-m
point-based pixel; CORINE is change area intersected with a 1-km
radius, and the epochs only overlap partially. NDVI green signal can
change without 411/412 mapped wetland loss due to wetland vegetation,
mowing, vegetation succession, weather, agricultural influence or
remote-sensing artifacts. It cannot identify which process occurs.


## Executed original data results (2026-11-20)

[Source-verified successful workflow #44673292597](https://github.com/zuizui0223/alan/actions/runs/44673292597), signed artifact #14600139114 (ZIP SHA256 ee11488f6eb0179acc8aacfc6ee028d5d6219c2ea5131a54e4e4f9721dd0c37b).
The exact 17-source-site output and 2500 block-bootstrap draws summary are committed at [17-site signed NDVI JSON](../results/all_species_forest/wet_falsification/signed_wetland_ndvi_17_20261120.json).

- Long-term signed OLS NDVI slopes were negative at **8/17** sites; robust Theil–Sen slopes negative at **10/17**; **8/17** had both negative.
- Direct early 2000–2004 versus late 2019–2023 full-growing-season medians fell at **13/17** sites. Median early-to-late change was **−0.01418 NDVI**; mean **−0.01654**.
- Mean signed OLS slope **−0.00409 NDVI/decade**, median **+0.00120**, indicating heterogeneity and dependence on summary. Spatial-block (8 independent blocks) bootstrap CI for the mean **[−0.01565, +0.00242]**.
- No sites had a mapped 411/412 CORINE wetland-class loss despite some changing NDVI; 4/17 had any 2000–2018 land-cover-class change.
- Exploratory Spearman between NDVI slope and candidate-minus-reference long-term *survey recording-rate* slope: **ρ=+0.147**, block bootstrap 95% **[−0.313,+0.657]**. This does not identify a coherent vegetation-mediated relative decline.
- Site ID 4297 shows the starkest example of differing spatial measures: CORINE within 1-km radius classified 60.52 ha from mixed forest 313 to peat bog 412 (2012–2018), while QA-screened 250-m-at-point NDVI OLS slope was **−0.0491/decade**, and early-to-late median changed **−0.1086**. Different footprints and epochs preclude interpreting this as verified bog degradation or proof against a physical expansion.

This does **not** support assigning a causal light effect; only one of 17 sites had long-term satellite brightening, the candidate plants were selected from partly the same discovery-era survey records, and 187-matched-pair geography remains unavailable. It establishes that mapped wetland-class transitions and vegetation greenness are not interchangeable measurements.
