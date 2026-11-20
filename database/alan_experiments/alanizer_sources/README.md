# ALANizer: recoverable field biomass and treatment-label audit

2026-11-20. No effect estimates or new pooled results. This audit recovers biological data linked by a methods report; it does not count the methods paper as an independent biological experiment.

## Primary sources and versions

- Read full bioRxiv v2 Methods/validation Results: https://doi.org/10.1101/2025.08.30.673291 ; source XML https://www.biorxiv.org/content/early/2026/02/19/2025.08.30.673291.source.xml . Hardware field deployment: 12 hedgerow sites at UBC Vancouver, at least50m apart, June–November2023; white/amber/dark cycling in5-day phases over20days. The two lit colours alternate with dark phases, with staggered site schedules. Light at1.5m produced approximately44/43lux below white/amber lamps and1lux at5m in calibration.
- Same-title final methods paper located through the author institution: https://doi.org/10.1111/jen.70181 . Final full text not read; bioRxiv API's published field remains NA, which is not evidence of no publication.
- Separate biological report located through author institution: https://doi.org/10.1016/j.biocon.2026.111798 . Abstract only read. It describes10,211 pitfall-captured arthropods and white/amber/dark treatment, reporting stronger capture under both light colours. That number is not the independent replicate count. Final Methods/Results have not been verified. Related SSRN5538162 PDF returned403; no bypass attempted.
- Original deposit: https://zenodo.org/records/15170657 . Metadata says June–December2023 and warns that planned labels can differ from actual treatment. Small files `metadata.csv` and `biomass.csv` downloaded and verified against repository MD5s. Large images/classifications and sensor series were not downloaded. Site/time identity supports a linked experimental programme; final-paper exact sample/exclusion correspondence remains pending.

## Verified data support

`audit.py` joins the two source tables by Site and StartDate, checking unique keys before merging. All134 biomass rows match metadata's actual Treatment. Eight of those rows disagree with the planned PrescribedTreatment. Thus using scheduled treatment would misclassify eight observed site-period records. The biomass table itself already carries the corrected treatment: this is not a correction to the author's numerical table.

There are12 sites;11 have all3 treatment arms in the biomass table. Site11 has no white-light biomass observation. The134 rows represent repeated site-periods with near/far measurements, not134 independent illumination units. Total equals Near+Far within floating-point precision. Start dates range2023-07-01 to2023-11-28. Units of biomass have not been verified and are not guessed.

The schedule metadata has396rows (12sites×33dates),249 with missing actual Treatment. Of262metadata rows without biomass,13 have a known actual treatment. Neither missing treatment nor missing sample is a measured zero. The records do not prove the reasons for exclusions or that every unsampled scheduled visit was attempted.

## Remaining requirements for numerical synthesis

Read the biological paper's final Methods/Results, establish mass units and trap-position mapping, verify exclusions and true exposure, and account for repeated dates and switching treatments within sites. Near/far traps share illumination; taxonomic levels would also share specimens. Do not count orders/families or near/far traps as independent experiments. A taxonomic predator category or capture biomass cannot establish predation rate. No wet/dry contrast or moisture moderator is inferred from the hedgerow setting.

Raw sources remain at `../gbif_annual/alan_experiment_sources/alanizer/`. Outputs retain the134-row crosswalk,8 planned/actual differences, site support and hashes. Reproduce with `python database/alan_experiments/alanizer_sources/audit.py`. Existing108-effect catalogue and three pooled results unchanged.

A second automated reviewer independently confirmed these counts, unique keys, treatment agreement and biomass sums from the original CSVs. This is not human double extraction. “Actual treatment” here means the author's corrected metadata, not independent validation against the unexamined illuminance time series. Carryover between treatment periods and capture effort also need checking before modelling.
