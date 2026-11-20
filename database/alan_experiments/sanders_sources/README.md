# Sanders2018 final-week field extraction

[Primary paper](https://doi.org/10.1016/j.cub.2018.05.078); [original EIDC data](https://doi.org/10.5285/d30168d3-6cbb-4d75-b73c-276e6083a1fe), Open Government Licence. Cite Sanders, Gaston and van Veen (2018) when reusing the dataset. Source URLs/hashes are recorded in `source_manifest.csv`.

The 2016 Penryn field experiment randomized eight light treatments within six blocks, with a separate light on each cage. Seven nonzero doses span 0.1–100 lux. Nine weekly surveys followed plant–aphid–parasitoid communities. Separate greenhouse and indoor experiments are not part of this extraction. Field plant-biomass records have lost labels and require a separate missingness audit.

## Calculated results and estimand

All 432 cage-weeks have unique keys; 48 cages remain at week9, one per block and treatment. Sampling factor equals one for every final-week record. Aphid totals include winged and unwinged columns; parasitoid totals sum the six host-associated mummy columns. Mummy counts are standing evidence of parasitism, not a rate per host or live adult abundance.

Fourteen final-week lnRR contrasts share six blocks and the same dark controls. Full covariance is retained (rank at most five). Approximate intervals use t5 with block delta-method variance. No sum of weekly observations is treated as extra replication.

Final-week aphid means range from24.17 at20lux to400 at50lux, versus56.67 in darkness. These differ from the authors' whole-season response and must not be presented as reproduction of that model. No monotonic dose response is imposed. This field programme is a new48-cage2016experiment, not the54 mesocosms used in Bennie2015/Anic2022.

Run `python database/alan_experiments/sanders_sources/extract_final.py`. Eligibility and second extraction adjudication remain pending, as for the other catalogue entries.
