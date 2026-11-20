# Soil isotope allocation: author summary data recovered

2026-11-21. Primary article: [Liu et al., DOI 10.1111/gcb.71029](https://doi.org/10.1111/gcb.71029). Author dataset: [Min Liu, ALAN dataset, version 1](https://doi.org/10.6084/m9.figshare.31288969.v1), CC BY 4.0. This acquisition supplements the existing [batch 6 primary review](../search/fulltext_batch6.md); it is not another paper or an additional independent experiment.

## Recovered material

The dataset–article association comes from the previously recorded primary Data Availability review in `search/fulltext_batch6.csv` (which explicitly lists this Figshare DOI), not from the Figshare API: its `references` and `related_materials` arrays are empty. The article-to-dataset link could not be independently rechecked against the full text on this pass and remains a provenance check to repeat when accessible.

The Figshare record contains one file, `Mydata_light_summary.csv`: 14 treatment–time rows and 17 paired mean/SE fields. Its downloaded byte size (3,296) and MD5 exactly match the API metadata. The unaltered file is retained here with attribution above. Time labels are 0, 1, 4, 8, 16, 32, 64, each with ALAN and Control. The file contains no plot IDs, sample sizes, units, data dictionary, covariance matrix, microbial network, or directly observed interaction records.

`extract_summary.py` preserves all 238 endpoint–arm–time records as text in `reported_summary.csv`, including two missing mean/SE pairs (`CUE13C` at time 0 in both arms) and two zero means (`Res13C` at time 0). These are not zeros added by extraction. `source_data_row` counts parsed data records, not physical text lines; the original CSV contains blank lines. No logarithms or pseudo-counts are applied.

## What can and cannot enter Results

The original screening identified five plots per arm and isotope allocation endpoints. This summary file alone does not verify the denominator of every SE, repeated sampling structure, units, or whether each field is isotope amount, recovery fraction, concentration, or a derived ratio. The variable names are therefore retained without expanding them into unsupported biological definitions. `Time=0` is not assumed to mean a pre-lighting baseline, and `Time=64` is not interpreted as an independent experiment.

In particular, `AGB13C` and `BGB13C` must not be treated as total aboveground and belowground biomass. A microbial correlation-network edge is not an observed feeding interaction. These data cannot establish stronger predation, reduced pollination, or a wet/dry habitat effect.

This acquisition adds **zero** effects to the current 108-row effect catalogue and **zero** records to the separately scoped 86-row reported-statistics registry. The local 238-row extraction is an additional source table, not 238 independent ecological effects. No effect variance is assigned from unverified sample sizes or zero covariance assumptions.

## Remaining verification

On this pass the PMC browser page returned a browser challenge and the official Europe PMC full-text XML route returned HTTP 500. Neither was bypassed. The earlier methods/results screening remains recorded; this is an access limitation on detailed rechecking, not a new exclusion or a claim that the article was never read.

Before quantitative reuse, recover the primary Methods/captions to map all units and time labels, establish SE denominators and sampling dependence, and obtain plot-level data or sufficient covariance for dependent contrasts. Check programme overlap with other Erguna lighting experiments. A summary-only author dataset is not proof that plot-level data do not exist.

Run `python database/alan_experiments/soil_isotope2026_sources/extract_summary.py` from the repository root to regenerate the long table and audit. The source hash is checked before extraction.
