# ALAN | Actual moth captures on islands versus historically prior artificial nighttime light

## Biological question, independently measured response

Can satellite nighttime light intensity explain **capture counts and community diversity of moths actually sampled on islands**, beyond photo-recorded plants whose species have only a positive GloBI moth-visitation literature claim?

This is an independent original moth-trap sampling frame, not iNaturalist botanical data and not evidence of effective pollen transfer. Positive phototaxis means **more moths in a light trap does not prove more living moths or more floral visitation**. Moth-trap method and local lamp placement are part of the detection process.

## Source verified 2026-11-10

Choi et al. 2026 Biodiversity Data Journal DOI 10.3897/BDJ.14.e211519, independent original Zenodo v2 DOI 10.5281/zenodo.21438943. The original archive has SHA/MD5-authenticated species-by-site counts and geographic metadata. Source source-quality results were computed from true original files in [ALAN GitHub Actions #42772784877](https://github.com/zuizui0223/alan/actions/runs/42772784877).

- 30 genuinely separate author-identified Korean coastal islands, 332 original trap location IDs (lat/lon), 377 unique island × site × survey-year rows and **30,442 actual captured moth individuals across 1,179 original species and 19 families**.
- Explicit true site-year survey dates are **17 years** drawn from 2007–2025 (no annual imputation). Original 2025 has **23 sites from only TWO islands**. Do NOT apply annual 2024 VIIRS to all 2007–2024 captures and call it contemporaneous or causal.
- Author's 2012–2025 subcohort has **323 genuine site-year observations on all 30 islands** (but sites/years are not independent). Only **44 locations** have measurements in >=2 different years; **19 islands** have >=2 surveyed years and only 12 have >=3.
- The original release provides a survey-year and capture assemblage but **no capture month, exact sampling date or per-site trap-night count**. Standardized method cannot erase unknown differences in actual trap effort, species attraction or season.
- Published species taxonomic source uses CP949 encoding. Unknown host plant is not an absent host plant and should never be imputed.

[Real survey count and temporal support ledger](../results/alan_real_korean_island_moth_sampling_support_20261110_verified.json).

## Historically correct night light exposure

For each original **island × site × true survey year**, retrieve the genuine original site pixel from *only the preceding calendar year's* Li & Cao v3 corrected PCNL (1992–2024, DOI 10.5281/zenodo.17013414). Actual survey years are 2007, 2008, 2009, 2010, 2011, 2012, 2014, 2015, 2016, 2017, 2018, 2020, 2021, 2022, 2023, 2024, 2025. Source rasters are 2006, 2007, 2008, 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2019, 2020, 2021, 2022, 2023, 2024 respectively. The 17 original, official-MD5-pinned raster years are split into four disjoint worker artifacts, then rejoined by exact original island+site+survey-year keys. Nodata/negative are missing, not zero-dark. This is a long-term annual PRE-survey exposure proxy, not a local streetlamp lux measurement.

The previous genuine 500m 2024 VIIRS test is available for the TWO 2025 survey islands as a geographically precise independent *pilot*. **Do not extrapolate** it retrospectively to all prior moth survey years. A 500m image from 2024 is a different product from cross-sensor calibrated PCNL, so coordinate matching alone does not make their values interchangeable.

## Statistical endpoints and exact hypotheses

Primary: **expected rarefied moth species richness at 20 captured individuals** per real site-year. Calculated analytically from original actual integer species counts using hypergeometric without replacement; drop sites with fewer than 20 actual moths instead of inventing richness. This removes simple dependence on the *number of captured individuals*, but does not establish equal trap-night effort or moth biology free of trap attraction.

Supporting endpoints: log(1 + captured moth counts), log(1 + actually observed species richness), rarefied richness at 10 individuals, and the Shannon exponential effective number of species.

Main exploratory contrast at same real island and same real sampling year: model endpoint on log(1 + PCNL_previous_calendar_year) and true source site elevation with **island × year fixed effects**, equal *total* weight per independent island, cluster-robust nominal uncertainty across the genuinely sampled island units. Require >=2 distinct geographic sites with different measured satellite values within a real island/year and >=0.02 range of original log1p previous-year PCNL. One island with many sites is not many islands. No synthetic zeros for unsurveyed years.

Hypothesis of interest:
- **H-aggregation:** Local light draws individuals into lamps or makes captures more likely; counts may increase while independent survival, effective pollen transfer or plant reproduction do not.
- **H-community filtering:** If light/catch method disproportionately favors some moth taxa, rarefied richness / evenness may respond differently from captured individual count.
- **H-landscape refugium:** A future true site- and year-matched dark-connectivity gradient may modify these recorded moth assemblages (the old 2024 VNL 5km maps cannot be assigned to older surveys without a time-consistent rationale).

All endpoints reported irrespective of coefficient direction or p-value. Sensitivities: site-year equal and island-year equal weighting, with/without site elevation. If actual within-island/year variation and island count are insufficient, mark the result NOT ESTIMABLE. Exact original source joins and declared 17 raster years precede any model.

## Interpretation guardrails

Light traps use positive phototaxis; trap count is a *sampling and behavior* outcome, NOT total moth population density. Pollen carriage and moth plant visits are unmeasured. No observational coefficient can prove light-caused biodiversity change, island-wide extinction, selfing evolution or effective night pollination. Historical moth field years have no precise dates/trap-night denominator, so apparent changes in captured composition can be seasonal/sampling. Additional independent site-year/date records or weather/land use are needed to argue mechanisms.

Status: genuine source original survey-year audit DONE; 17-source-year matched PCNL acquisition and original field-data fits are implemented and evaluated only after all original raster MD5/keys pass. Preserve original ALAN ROI9 coefficients and PR #34 geospatial result unchanged.
