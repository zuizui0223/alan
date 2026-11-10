# ALAN | Original 2012–2024 annual satellite dark-patch exposure vs true worldwide 4-family GBIF photo record composition

Design frozen 2026-11-10. This is an *observational record-composition* study, not a field moth census or a pollen-deposition experiment. Treat the source null and source-overlap failures as interpretable results.

## Exact source, geography and ecological estimand

- GBIF 2026-08 DWCA preview DOI 10.15468/dl.8yrbe7, source-verified all-records extraction #42782146115, original 22,713,063 georeferenced human observations 2013–2025 from original Noctuidae, Geometridae, Nymphalidae, Pieridae. 1,217 source-data GSHHG islands and 228,361 island×0.05-degree geographic photo-cell×year rows.
- Historic corrected PCNL v3 2012–2024 source DOI 10.5281/zenodo.17013414 assigned only to the original photo-OBSERVATION year's previous year, at the 0.05-degree original cell center **only if the center falls on the same named GSHHG island**. Preyear source SHA256 64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852, complete source #42786502922.
- The registered main original recorded-guild-comparison cohort: 196 independent islands, 1,083 observed island×year groups, 62,688 original photo-cell×year rows (both recorded guilds, >=5 original records, >=2 source light levels, >=0.02 log1p prior-year PCNL within island×year). Main result is photographic Noctuidae+Geometridae record share vs the same two moth and Nymphalidae+Pieridae butterfly family records, not actual natural moth abundance or flower visits.
- Original primary source model #42787302661: original independent-island-equal β light=−0.8651 percentage points per log1p prior PCNL [−1.9737,+0.2435] p=.1254; no resolved association. Source volume strongly correlated with record proportions and site-equal weighting β=−1.4754 [−1.9656,−0.9852], a DIFFERENT estimand driven by very uneven sampling; do not assert globally robust moth decline from that. Exploratory stage only.

## Original true annual dark habitat proxy

At each of those exact 62,688 source photo-grid centers, sample the **actual year−1 corrected original PCNL raster** in a geodesic 5km radius and the original GSHHG v2.3.7 high-resolution L1 island land polygon **minus original L2 inland lake polygons**. Never call sea or lake black pixel dark habitat. At least 16 valid pixel land-centers and >=80% raster coverage required.

Two source-prespecified PCNL-native darkness thresholds: zero and three. These are native intercalibrated satellite values, not calibrated lux at flower height or objectively verified nocturnal flight corridors. Zero detected-dark-pixel => largest-dark-component conditional fraction undefined, never zero. Source geometry and antimeridian problems explicitly unsupported. Threshold repeats on the SAME photo cell are NOT independent sample replicates.

Three geometry variables per site: local mean radiance in the 5km landmask; dark-land fraction; largest 4-neighbour-connected dark patch / all dark land. This is a 1km coarse global-satellite pattern, not individual 0–100m streetlight phototaxis. It cannot recover narrow real habitats or dark corridors under canopy.

## Conditional ecological SOURCE test

Main model: original documented moth-family GBIF photo-record percentage ~ yearminus1 PCNL pixel intensity + 5km dark land fraction + conditional largest connected dark patch fraction + log1p(all four-family source observation records) + geographic cell-center latitude, with **actual island×2013–2025 observation-year fixed intercepts, independent island-cluster covariance**, source-independent-island equal total weight. Select island×year groups with >=2 true photosites, >=0.02 local preyear loglight span and >=0.05 observed conditional dark-connectivity span. Recalculate support after excluding no-dark, missing and raster-unsupported locations; do not claim all original 196 islands remain identified.

Compare (1) original intensity-only source model on the same supported rows, (2) joint intensity+dark-fraction+dark-patch connectivity, (3) photo-documentation-volume as a negative-control outcome, (4) photo-site-equal weighting, (5) requirement >=0.05 dark-fraction contrast. If fewer than 30 independent source islands, design matrix lacks rank/overlaps with original PCNL enough to be singular, or condition number exceeds 1e11, report **not estimable**; no significance claim.

Source truth is GBIF human-photo documentation of 2 moth families and 2 butterfly families; lack of a GBIF photo is not verified absence. Moth species could be active by day; butterfly visitors are not a controlled experimental negative group. Year-, region-, photographer-/dataset-, urban landcover-, trap effort-, month-, and site-specific selection bias remain. A significant conditional correlation is NOT proof of light-caused pollination or moth abundance. A null could reflect coarse PCNL pixels, sparse true dark-patch contrasts, or observer bias.

At this stage only the source table is prepared. Full past-year 5km raster scan is independently executed under source run https://github.com/zuizui0223/alan/actions/runs/42788295115; do not declare model effect before the original source workers and joint analysis pass.