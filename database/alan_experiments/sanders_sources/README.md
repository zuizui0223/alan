# Sanders2018 final-week field extraction

[Primary paper](https://doi.org/10.1016/j.cub.2018.05.078); [original EIDC data](https://doi.org/10.5285/d30168d3-6cbb-4d75-b73c-276e6083a1fe), Open Government Licence. Cite Sanders, Gaston and van Veen (2018) when reusing the dataset. Source URLs/hashes are recorded in `source_manifest.csv`.

The 2016 Penryn field experiment randomized eight light treatments within six blocks, with a separate light on each cage. Seven nonzero doses span 0.1–100 lux. Nine weekly surveys followed plant–aphid–parasitoid communities. Separate greenhouse and indoor experiments are not part of this extraction. Field plant-biomass records have lost labels and require a separate missingness audit.

## Calculated results and estimand

All 432 cage-weeks have unique keys; 48 cages remain at week9, one per block and treatment. Sampling factor equals one for every final-week record. Aphid totals include winged and unwinged columns; parasitoid totals sum the six host-associated mummy columns. Mummy counts are standing evidence of parasitism, not a rate per host or live adult abundance.

Fourteen final-week lnRR contrasts share six blocks and the same dark controls. Full covariance is retained (rank at most five). Approximate intervals use t5 with block delta-method variance. No sum of weekly observations is treated as extra replication.

Final-week aphid means range from24.17 at20lux to400 at50lux, versus56.67 in darkness. These differ from the authors' whole-season response and must not be presented as reproduction of that model. No monotonic dose response is imposed. This field programme is a new48-cage2016experiment, not the54 mesocosms used in Bennie2015/Anic2022.

Run `python database/alan_experiments/sanders_sources/extract_final.py`. Eligibility and second extraction adjudication remain pending, as for the other catalogue entries.

## Field biomass allocation and missingness audit (2026-11-20)

The original bean and barley files each contain30 of48 cage records;23 cages have both outcomes and11 have neither. The author RTF describes harvest at experiment end, one pot per mesocosm, and missing labels. Missing values are not mortality or zero biomass; random missingness is not established.

**Cage39, block5 conflicts across files:** the population table consistently assigns5lux, whereas both biomass files assign1lux. All retained biomass block IDs agree with the population table, but this light label does not. Neither file is silently treated as the corrected allocation. The existing14 insect contrasts remain conditional on population-file treatment labels; their cross-endpoint allocation validation and final adjudication are now explicitly unresolved.

Using biomass-file labels only, complete lit/dark block support ranges0–3 for bean and2–5 for barley. Bean10lux has no complete dark-paired block. These counts describe available data, not an eligibility rule or a selected effect analysis. No biomass effects are added, no missing values are imputed, and insect effects are not relabelled or recomputed. A corrected author allocation map is needed before joint insect–plant interpretation; label loss must also be addressed before inferring a full-experiment biomass response.

Reproduce with `python database/alan_experiments/sanders_sources/audit_biomass_support.py`. Outputs preserve all48 cage IDs per species, the two source labels, arm availability, complete-pair support, conflicts and archive/member hashes. Source: [original author archive](https://doi.org/10.5285/d30168d3-6cbb-4d75-b73c-276e6083a1fe), `ALANintensityDescription.rtf` and three field CSVs. Exact coordinates remain unverified; the archive identifies Penryn campus but this is not a measured cage location.

A separate agent independently reread the ZIP and confirmed the48-cage design, both30-record biomass files,23 shared and11 jointly missing cages, the Cage39 conflict, and all paired-support counts. Zero complete pairs for bean10lux does not mean there are no10lux biomass records. This is an AI computational/source check, not human double extraction.
