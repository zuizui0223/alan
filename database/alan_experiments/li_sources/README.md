# Li2026: reproductive outcomes, held outside the main catalogue

Primary preprint: https://doi.org/10.64898/2026.07.10.737511
Public author data: https://github.com/Dongboli/experimental-data/tree/main/Artificial%20lights%20at%20night

The public README labels the data for peer review only. The original table is cached locally and is not redistributed here. `audit.json` records its Git blob and SHA256; the extraction script verifies the blob before calculation.

The paper describes three field blocks, dark, point and array lighting, with two endpoint plants per treatment and round. Lamp number also differs between lighting arrangements. Pollinator visits and day/night contributions were not measured directly. This is a preprint.

## Independent data audit

The deposited table contains 79 rows versus 90 planned plants. ML and ZC have unresolved counts or label discrepancies; one seed outcome is missing. Missing outcomes are not zeros. Only SN (mapped to Silene nutans from the listed species initials) passes the complete-cell rule: 18 plants, two per treatment in each of three blocks. Selection uses design completeness, not effect direction.

`SN_descriptive_effects.csv` contains six correlated descriptive contrasts, not six independent experiments. Endpoint plants are averaged within each block and treatment. Means, paired block covariance, delta-method variances and approximate t2 intervals are retained. Probability intervals are unconstrained Wald intervals and can extend beyond feasible bounds; they are not bounded probability estimates. Three blocks provide weak precision.

For seeds, array/dark and point/dark mean ratios are 6.903 and 2.030, respectively; both approximate lnRR intervals include zero. These positive point estimates are retained as counterevidence, not evidence of established improvement or moth-mediated effects. The six rows remain outside the 66-effect main catalogue pending whole-study adjudication.

Run `python database/alan_experiments/li_sources/extract_complete_species.py` with the pinned original table in the external source cache. Outputs include cell counts, block means, six effects, their joint covariance and audit.
