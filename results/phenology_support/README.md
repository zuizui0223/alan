# Flowering phenology and sampling overlap in the frozen ALAN frame

The comparison uses the primary analysis's 78 five-degree regions, 10,698 site-years and 195 region-years in 2019, 2021, 2023, 2024 and 2025. The 544 shared WCVP plants occur in 72 of those regions. The primary island-minus-mainland estimate of −1.902 and the secondary 64-species pollen-limitation analysis keep their existing status.

## Sampling opportunity

Within each region-year, classify sites above the median site-year VIIRS radiance as bright and the remainder as dark. The same threshold applies on islands and the mainland. These are relative light environments. Site-year radiance is not the light level at an individual flower.

Count distinct local observation dates for the same species, region and year in each of the four land-class and light combinations. Then repeat within calendar month. As a background effort proxy, count distinct dates with any plant photograph in the full frozen corpus (61,096 photo records) and require two or more such dates within each of the four arms during that month. This measures calendar coverage. Observer-hours, visits without photographs and systematic searches for the focal species are unavailable.

- All four arms in the same year contain 65 photos of 10 species in three regions. Nine species have just one distinct date in the least-sampled arm. The tenth has two.
- All four arms in the same month contain three species in two regions. All three pass the background two-day proxy. Each has only one observation day of the focal species in its least-sampled arm that month. No species has a second common month within that region-year.
- Background effort alone meets the two-day rule in at least one month in 26 region-years in 13 regions. It meets that rule in three or more months in nine region-years. These background days provide no presence or nonflowering records for the focal plant.

## Evidence for flowering and dates

The public iNaturalist API returned observations for 8,457 of 8,466 original lookup identifiers (8,444 distinct public observation UUIDs). The frozen photo ID and local observation date both matched for 8,456 records. Of these photo-linked records, 5,547 had an observation-level Flowering annotation and 69 had No Evidence of Flowering. One observation can carry multiple phenology terms. Neither an unannotated observation nor a different phenology term supplies a nonflowering observation. Counts here are photo-linked records, not independent individuals.

All 65 focal images from the four-arm region-years were visually scored. Open flowers of the focal plant were visible on 39, none were visible on 15, and 11 were unresolved. The second category describes what the image shows and cannot establish a whole plant's flowering status. The public original-image headers were checked for EXIF DateTimeOriginal (first 256 KiB). None of the 65 yielded the tag. Dates used here are observation dates. The manual review is restricted to all four-arm candidates. Other images have only the public metadata check.

No species-region-year had even three distinct sampling dates per arm. There is no time series with comparable coverage across the four arms for resolving a local flowering peak and both ends of flowering. We therefore report no island–mainland contrast for peak shift or duration. First and last uploaded or observed flowers would confound phenology with sampling.

## Reproduction

Run scripts/analysis/alan_phenology_frame.py with the frozen primary inputs directory, signed inputs directory and a private output directory. Run alan_phenology_annotations.py and then alan_phenology_support.py on that same output directory. The former caches public API responses locally. The latter matches manual review labels by SHA-256 of the frozen photo ID, then writes the aggregate results beside this file. Public annotations can change with time. The archived tables describe the 2026-11-12 check.
