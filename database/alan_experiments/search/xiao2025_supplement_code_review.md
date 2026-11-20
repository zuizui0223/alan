# Xiao2025 supplement and code audit

Reviewed 2026-11-20. No effect extracted; numerical extraction remains pending.

## Primary-source summary (under 200 words)

[Primary Methods/Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC12152749/) describe three paired light/dark compartments for plant growth, with a separate fourth pair for feeding assays. The greenhouse coordinates are 125°24′03″ E, 43°59′49″ N (125.400833, 43.996944). Pots are subsamples, not independent lighting replicates. Final biomass was collected from 284 surviving pots on 7 December 2023. The feeding assay has only one independently treated compartment per arm.

The six-page [supplement](https://doi.org/10.6084/m9.figshare.29118899) contains two variance-structure tables and three figures. Tables S1–S2 report model standard deviations/multipliers, not treatment-arm sampling dispersions. Figure S1 supplies photographs and a factorial schematic, not numbered compartment assignments or pair IDs. Figure S2 biomass panels summarize herbivory/nutrient effects; its ALAN panel concerns root mass fraction. Figure S3 presents photosynthesis means and SEs, including negative means. No compartment-level mean, sample-size or covariance table is supplied. The supplement title-page DOI contains 2025; the verified article DOI is [10.1098/rspb.2024.2811](https://doi.org/10.1098/rspb.2024.2811).

## Independent code and extraction audit

Inspected [Zenodo R code](https://doi.org/10.5281/zenodo.14000092), cached as `R_code.txt`, alongside all supplement pages, including rendered figures. The code loads `ALAN_data.xlsx` and `relative_growth_rate_data.xlsx`; neither workbook is embedded. Dryad metadata lists these two files and README under [10.5061/dryad.1rn8pk143](https://doi.org/10.5061/dryad.1rn8pk143). No raw observations or executable inline data replace them.

- Plant model code names `Compartment`; larval code names `Region`. No compartment-to-pair lookup is defined. Geographic coordinates cannot recover this allocation. Obtain workbook and README before deciding paired covariance or aggregating pots.
- `summarySE` calls group by treatment/status (and sometimes nutrient/herbivory), without first averaging within compartment. Their SE bars must not be treated as independent-light-unit dispersions, even if digitized. The plotted summaries cannot reconstruct compartment covariance or pair differences.
- S2 biomass means pool lighting levels. S3 means are potentially digitizable descriptive information, but negative photosynthesis makes ordinary lnRR invalid, independently of the replication issue. No digitization was performed.
- `mydata$Photo` and `mydata$Tolerance2` are referenced although no `mydata` assignment appears. A clean-session rerun requires checking this apparent object-name error against raw inputs; no correction or rerun was attempted.
- The separate feeding assay must not gain independent-light replication from larvae, plant species, or four subcompartments.

## Species crosswalk check (own database audit)

All 12 binomials and China native/invasive statuses in `xiao_sources/species_as_reported.csv` match primary Table 1. Independently rereading EIV and the wet-comparison table reproduced four exact-name numeric MOIST values: Amaranthus retroflexus 4.6, Artemisia annua 4.1, Ambrosia artemisiifolia 4.0, Cenchrus spinifex 1.5. Lolium perenne and Datura stramonium have `x`, not numeric dry affinity. The other six names lack exact EIV matches. Intersections with both wet-comparison source names and accepted names are empty (363 labels, 361 accepted species). This exact-name check is not synonym resolution and cannot establish ecological absence or infer treatment response.

Verdict: greenhouse growth design remains eligible in principle, but neither recovered code nor supplement supports independently replicated effect-size/variance extraction. Do not convert model variance multipliers, pot SEs, or reported significance into compartment-level effect variances. Keep feeding-only estimates outside replicated-light pooling.
