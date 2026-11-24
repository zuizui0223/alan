# Water quantity, proximity and patch-geometry source gate (2026-11-24)

## Executed exploratory macro reanalysis
Fixed 1706 valid sites from frozen 1707 coordinate panel and baseline 1999 GLAD distance; source blobs in JSON. Original area0 is the first two available survey-year means of percent pixels with water detection >=50% in approximately 1 km square, not the same exact 1999 disk or an actual wetland extent. Frequency, area0, changes, weather, built cover, observation, group richness controls retained. Distance N=1-distance_m/1000; water W=log2(1+area0); L=light/10. L,N,W,LN,LW,NW included, and separate model adds LNW. Same block demean and block-cluster score variance as original, with t quantile approximation. FWL re-orthogonalization; 4 overlapping scopes (not independent tests): all, detected within 1km, detected and area0>=1%, latter with positive light change. Report all outcomes/models; nominal CIs, not multiplicity-corrected. Conditional slopes describe recording-rate trends, not abundance or causal effects.

Common support: near <=100m/high area>=1% 143 sites/42 brightened; near/low 5/0; far/high 335/85; far/low 1223/204. Pearson correlation N/W .8115 overall, .7378 detected. No four-cell balanced brightening contrast. The 1% classification comes from prior wet_water source; 100m from prior distance analysis, not fitted breakpoints.

| Scope | Used sites | LW candidate minus reference | 95% approximate CI | LNW candidate minus reference |
|---|---:|---:|---|---:|
| all | 1674 | -0.779 | -2.083, 0.526 | 1.357 |
| water_detected | 846 | -0.561 | -1.962, 0.840 | 3.216 |
| water_rich | 458 | -1.043 | -3.478, 1.392 | 4.659 |
| bright_water_rich | 117 | -2.445 | -8.501, 3.611 | -6.343 |
Water-rich conditional LW slopes both negative (candidate -3.640 [-6.413,-.867], reference -2.597 [-3.757,-1.436]); differential unresolved. All LNW intervals include zero. Bright-rich effective information as low as 2.56 blocks for LNW; wide intervals. No support for a candidate-specific distance/quantity amplification. No true connectivity estimate or new experimental/meta coefficient.

## Original geometry source gate
Three prior coordinate-bearing source tables, exact ID validation (tolerance 1e-7 degrees), 109 unique IDs overall, only 59/1706 in current valid panel, 6 brightened, 15 blocks, 2 near100m. No original annual.duckdb or raw tiles present in current session; direct network DNS unavailable. File search returned no matching local source; not proof source does not exist elsewhere. Do not reconstruct coordinates from block or distance. Do not treat these selected 59 as representative. Full local extraction requires original geo table and cached 1999 GLAD tiles.

New local extractor validates all 1706 coordinate IDs and reproduces each original nearest-water distance (1e-3 m tolerance), then extracts 250/500/1000 m disks, annual detection thresholds 50/80, four- and eight-neighbour component size metrics. Water quantity and same-patch probability sum(p_i^2) distinct. No water -> undefined topology; unknown pixels -> topology undefined (require complete disk); water fraction withheld below 95% coverage. Disk clipping and pixel resolution affect components. These measure open-water patch geometry only, not stream-network continuity, wet vegetation, or species-specific insect movement. Original coordinates never written. Core patch_metrics function and synthetic tests were executed successfully in local Python: identical 8/81 pixel water fraction and nearest distance one pixel, connected probability 1 vs eight separate patches .125; diagonal adjacency, missingness and no-water cases passed. Full real-raster extractor and repository standalone Node replay not run in Actions. No real patch metrics estimated.

## Reproduction
node scripts/analysis/alan_water_quantity_distance_20261124.mjs
python scripts/analysis/alan_water_patch_geometry_20261124.py --self-test
Full extraction command in extractor docstring; requires original local inputs. No main changes. No claim of identified causal plant mediation.
