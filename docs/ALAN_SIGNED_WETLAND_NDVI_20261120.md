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
