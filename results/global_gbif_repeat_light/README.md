# Worldwide GBIF repeat-coordinate screen

The 1 November 2026 worldwide GBIF release contained 3,887,364,403 occurrence rows, including 621,820,608 vascular plant rows. Every one of the 10,346 source partitions was scanned. The private project database links 67,105,132 eligible records to their source dataset, accepted species key, basis, recorded coordinates and date. The snapshot is [GBIF download 10.15468/dl.jeesns](https://doi.org/10.15468/dl.jeesns). The ten sources contributing to the final survey comparison and their citations are in [survey_sources.csv](survey_sources.csv).

## Cohort and interpretation

The two windows are 1992–1999 and 2020–2024. A repeat requires exactly matching published coordinates, source dataset and basis in both windows, and the species itself in both for the persistence screen. The separation is at least 21 calendar years. Endpoints do not require 20 annual observations and do not exhaust all possible 20-year windows. Rounded coordinates and reported uncertainty up to 1 km can still encompass different individuals or observer routes. Complete calendar fields, PRESENT records and human observations, observations or preserved specimens were required. Jan 1 was excluded from the survey reconstruction. Sites where the species disappeared entirely are outside the contrast.

There were 18,588 repeated source locations, 156,894 repeat species–source-location combinations and 2,966 species. A proxy survey day requires five distinct vascular plants at the source location. The same calendar month needs two such days in each epoch. This retained 351 survey days, 64 location–month combinations, 2,041 species–location pairs, 650 species, and 58 distinct coordinates. Of these coordinates 15 were on two islands, 42 on the Eurasian mainland and one on the North American mainland. The source scan was worldwide. The final comparisons represent a heavily European sample. See [coverage](coverage.json), [raw repeat species](repeated_taxa.csv), [geography](geography.csv) and [all species](species_screen.csv).

The outcome is the percentage of eligible survey days containing a species. Period rates are calculated by month with equal month weights. Satellite brightness in each period is matched to the actual years of the qualifying visits, averaged in the same monthly strata. Each physical coordinate has equal weight in its species' regression. The model compares changes within 100-km equal-area blocks controlling for initial brightness and change in mean list richness. Dataset-date lists proxy recording effort. Observer turnover, list targeting, habitat change and survey duration remain unresolved. A brighter satellite pixel can represent lighting outside the plant population. [Li et al. harmonized DN, v10](https://doi.org/10.6084/m9.figshare.9828827.v10)

## Screening outcome

Formal screening requires at least 20 coordinates, ten blocks with information on the adjusted light gradient, and adequate brightened-site replication. The two negative contrasts are the gradient of recording change on brightening and the mean recording change on brightened sites. No species supported both tests. No multiple-testing-adjusted decline or island interaction could be estimated.

Only two taxa met the initial 20-coordinate and ten-block counts. Descriptive coefficients below use the same adjustment with the available blocks. They are unstable and have no inferential status. [Complete descriptive comparison](high_coverage_descriptive.csv)

|Species|Island / mainland|Mean change at brightened sites (percentage points)|Mean change elsewhere|Adjusted slope per +10 DN|Informative blocks|
|---|---:|---:|---:|---:|---:|
|Plantago lanceolata|7 / 13|−31.88|+6.32|−34.46|5|
|Dactylis glomerata|8 / 13|−32.92|−4.17|+40.67|3|

## Flowering follow-up

Plantago lanceolata was the sole directional lead among the two best-covered species. The Chronicles of Nature Calendar contained 21 quality-accepted flowering onsets from three mainland monitoring locations in 2004–2014. All 21 had a same-year satellite estimate. Flowering followed landscape green-up by 9–61 days. Harmonized light was 0 DN in all 21 observations. No location matched a GBIF candidate coordinate exactly. There is no brightness contrast in this panel and no island flowering series in it. The daylight or pollinator consequences of these intervals were not measured. The 23-species European PEP725 herb and dwarf-shrub release also screened had no P. lanceolata. The UK first-flowering paper points to Nature's Calendar for the observation data; its supplement alone does not provide a matched island panel. The coverage verdict applies to the sources examined.

[Field coverage](exploratory_field_coverage.csv) · [Chronicles of Nature Calendar v1.0.5](https://doi.org/10.5281/zenodo.3607556) · [NDVI phenology and method](https://doi.org/10.1038/s41467-019-10235-8) · [PEP725](https://pep725.eu/dataset/) · [PEP725 onset extract](https://zenodo.org/records/21979471) · [UK flowering study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8808087/)

## Replay and scope

Run alan_global_gbif_source.py on the dated snapshot, then alan_global_gbif_database.py, alan_global_gbif_exposure.py, alan_global_gbif_screen.py and alan_global_gbif_field.py. These scripts are under scripts/analysis. The DuckDB file and the adjacent snapshot Parquet shards are held together locally under work/gbif_world. Moving the source folder requires refreshing the occurrences view. Field data are actual onset dates matched to GIMMS3g spring green-up (method 4) and harmonized night light at the field coordinates in each year. The two older manuscript results retain their roles, with −1.902 primary and the 64-species pollen-limitation association secondary.
