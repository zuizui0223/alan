# Swiss plant observations and satellite night light

## Scientific questions

1. Which plant species changed in detection between 1992–1999 and 2020–2024 in comparable GBIF botanical lists? Did those changes follow local brightening?
2. Within a phenology station, do brighter years have a larger absolute separation between floral onset and satellite spring green-up? The primary phase is the reported start of flowering; the separate 50%-flowering analysis extends further back and includes herbs. The signed interval (flowering day minus green-up day) and its components retain timing and direction.

## Occurrence data

[Info Flora / InfoSpecies, Swiss National Databank of Vascular Plants](https://doi.org/10.15468/7jffp3), CC BY 4.0. Source Darwin Core archive for the GBIF registered dataset UUID 83fdfd3d-3a25-4705-9fbe-3db1d1892b13, release dated 2026-10-22, 16,302,425 rows. Human observations in Switzerland only; present, complete dates, vascular plants, with positive or zero reported uncertainty up to 3,536 m. The provider generalizes most records to 5-km squares. Species are formed from genus and specific epithet, rolling up subspecies/varieties. Coeloglossum viride joins Dactylorhiza viridis, and Anemonoides nemorosa joins Anemone nemorosa by exact GBIF synonym. Other spellings are retained. Records explicitly citing iNaturalist in the selected provenance fields are excluded; undocumented overlaps with other sources remain possible.

A list is a proxy comprising the distinct plant species recorded by the same named recorder on one date in one Swiss national 5-km square. Each list must contain at least ten species. Duplicate species within a list collapse; blank or generic observer labels and January 1 dates are excluded from the list comparison. Restrict to the same recorder string, square and calendar month in both periods. For each species compute the change in percentage of lists containing it in each matched stratum and average across strata within the square, including matched strata where the species was unrecorded. Across cells that recorded the species at least once in either period, estimate the equal-cell mean change and the light association. The light coefficient compares cells within the same 25-km block, controlling initial light and the change in list richness. Light is the 25-point average within each 5-km cell in each survey year, averaged over matched survey dates. At least 20 recording cells and ten spatial blocks are required. SEs cluster on 25-km blocks, with t reference; BH correction across species separately for each outcome. Survey duration, targeted visits, recorder identity changes behind a single name, and within-cell location changes cannot be reconstructed completely. List detection measures neither population size nor confirmed local extinction. Two endpoint periods cannot establish uninterrupted annual decline.

## Flowering and NDVI

[MeteoSwiss open phenological observations](https://opendatadocs.meteoswiss.ch/a-data-groundbased/a9-phenological-observations) provide actual calendar dates for station-monitored flowering. All 175 published station files were checked. Match station coordinates to annual GIMMS3g NDVI spring onset, 1/12°, derived as passage through 50% of that pixel's fitted annual NDVI amplitude (Method 4), and to the same year's harmonized satellite night light. The primary signed interval is measured in days, positive when flowering follows green-up. The absolute interval tests widening regardless of event order. These capture the relation between a monitored species' flowering and surrounding vegetation, not measured harm or loss of overlap with a pollinator. The station coordinates describe the monitoring area, and satellite grid cells describe a landscape. [Wang et al. (2019), GIMMS phenology](https://doi.org/10.1038/s41467-019-10235-8) · [Data and fair-use description](https://data.globalecology.unh.edu/data/GIMMS_NDVI3g_Phenology/Readme_and_Fair_Data_Use_Policy_GIMMS-Phenology.pdf) · [Harmonized night light, Li et al. v10](https://doi.org/10.6084/m9.figshare.9828827.v10).

A station–species combination needs at least three matched years. Estimate each species with station and year fixed effects, SEs clustered on Swiss national 25-km blocks and finite-cluster t intervals. At least ten stations and ten blocks per species. BH across species per outcome and phase. The pooled estimate is an equal-species-weighted common slope with species×station and species×year fixed effects and the same spatial clustering. Coefficients refer to 10 harmonized DN units, not lux. Exposure can saturate, and local climate and land-use change can covary with brightness. The primary analysis uses all eligible field observations. A second result requires an independently sourced GBIF occurrence of the same named species in the station's 5-km grid and year; this supports local occurrence but does not identify the monitored individual. No GBIF occurrence date is used as a flowering date. Source coverage is for Switzerland, not the island–mainland contrast. Historical −1.902 and the secondary 64-species pollen-limitation result retain their existing interpretation.

## Outputs and replay

- [GBIF coverage](gbif_coverage.csv), [species comparisons](gbif_species.csv), [MOIS crosswalk and trait evidence](gbif_mois_evidence.csv), [three MOIS candidates](previous_three.csv).
- [Coverage by plant, phase and overlap](phase_coverage.csv), [equal-species timing estimates](pooled_intervals.csv), [onset](flowering_models.csv), [50% flowering](full_flower_models.csv), and [GBIF-confirmed onset](gbif_overlap_onset.csv) and [50% flowering](gbif_overlap_50pct.csv).
- Raw records, observer identifiers and locations are not redistributed. The source checksum is recorded separately. The GIMMS provider requests early contact when preparing a publication that uses these data.

Download the source archive via the GBIF dataset's DWC_ARCHIVE endpoint. LIGHT contains the previously restored annual harmonized DN GeoTIFFs, using calibrated DMSP through 2013 and simulated VIIRS after 2013. From scripts/analysis run:

    python alan_gbif_swiss_records.py CACHE/swiss.zip CACHE
    python alan_swiss_phenology_sources.py CACHE/pheno
    python alan_swiss_flower_gap.py CACHE/pheno LIGHT OUT
    python alan_gbif_swiss_lists.py CACHE LIGHT OUT
    python alan_gbif_swiss_prior.py OUT
    python alan_swiss_gap_summary.py CACHE OUT
## Findings from this release

Of 10,867,100 eligible Swiss plant occurrences in 1992–2024, the published lifeStage field was empty for all. The 60,385 qualifying lists in the two endpoint periods yielded 350 comparable lists (180 early, 170 recent), 79 recorder × square × month pairs, 28 recorder labels and 59 squares. Seventy-two species met the temporal comparison requirements. Achillea millefolium declined in matched list detection by 52.1 percentage points (q=0.0483). Its light comparison had only three spatial blocks with residual light variation after covariate adjustment. For the light comparison we count blocks actually contributing residual light variation and require at least ten. Eleven species met that criterion across the screen; none had an FDR-supported negative light association. The decline of this single species cannot identify shared reproductive traits. Of the previous three MOIS species, two orchids appeared in the matched lists but with too few squares and the third had no eligible source record. Absence from these data does not establish biological absence.

Equal-species estimates for widening of the absolute gap, in days per 10 DN (95% spatial-cluster intervals):

| Field phase | All eligible stations | Same species independently recorded by GBIF in the grid-year |
|---|---:|---:|
| Start of flowering | −0.43 [−1.99, 1.14], 13 species / 86 stations | +2.44 [−5.38, 10.26], 4 species / 34 stations |
| 50% flowering | +0.06 [−1.34, 1.47], 22 species / 87 stations | +0.87 [−3.88, 5.62], 8 species / 65 stations |

The complete onset panel spans 1996–2014, the complete 50%-flowering panel 1992–2014. The GBIF-confirmed panels are selected subsets of these station observations, and their wider intervals permit sizable effects in either direction. No species in the full station onset or 50%-flowering panel had an FDR-supported widening. The pooled estimates use species × station and species × year fixed effects. A separate explicit-dummy regression reproduced both Prunus avium single-species absolute-gap estimates and clustered intervals. Species and stations contributing to the pooled fit are counted after minimum repeat-observation and spatial support checks. Coordinate generalization, observational effort, and residual weather/habitat confounding limit interpretation. These data do not measure moth activity, effective pollination, or seed production.
