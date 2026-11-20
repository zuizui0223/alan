# Zeale2018: paired site feeding-attempt contrasts

[Author data](https://data.bris.ac.uk/data/dataset/31dvq1elivhby2dap2tm0zgz94), DOI10.5523/bris.31dvq1elivhby2dap2tm0zgz94; [primary design review](../search/zeale2018_design_review.md). Retrieved2026-11-20 directly from the data server after the DOI redirect failed. The archive supplies activity, telemetry and a short README. source_manifest.json records URLs, downloaded sizes and locally computed SHA256 fingerprints (not independently published checksums). Raw files remain in the external cache.

## Extraction

The320-row activity table contains eight sites, five conditions, two hedge sides and four bat groups, with unique keys. Total passes57558 matches the paper. Feeding data occupy only40 rows: Pipistrellus on the designated lit side. Missing feeding observations elsewhere are not zero. Recorded ratios are rounded; reconstruct each site-night ratio from its observed buzz and pass counts, all with positive denominators. The320-row telemetry file was acquired but not used here.

For each endpoint and light colour, compute log(mean of eight light-site values / mean of eight corresponding dark-site values). Variance is var(light/mean(light) - dark/mean(dark))/8; intervals use the delta approximation and t7. For buzz/pass, the endpoint is the unweighted mean of site-level ratios, NOT the pooled ratio of total buzzes to total passes. No colour was selected according to its result. All12 contrasts and their full12x12 covariance are retained, including shared dark controls, endpoints and sites. Matrix rank is at most7; do not invert it or treat12 contrasts as independent studies.

| Colour | Pass ratio | Buzz ratio | Mean site buzz/pass ratio [approx95% interval] |
|---|---:|---:|---:|
| Orange |8.574|8.113|1.145 [0.791,1.657]|
| White |5.521|5.454|1.128 [0.519,2.455]|
| Green |11.246|10.701|1.237 [0.643,2.379]|
| Red |2.472|2.371|0.978 [0.747,1.280]|

These are descriptive paired mean contrasts, not a reproduction of the original negative-binomial/logistic models. Our individual unadjusted intervals for red pass/buzz counts exclude1, whereas the original post-hoc model comparison did not identify a clear red effect. Different estimands, models and multiplicity handling mean these should not be substituted for one another. No colour-specific threshold claim is made. The ratio intervals including1 do not establish equivalence or lack of biological effect.

## Limits and adoption

Darkness was always measured first. Pairing controls stable site differences but cannot separate lighting from subsequent night/order changes. Opposite sides are not independent control sites. Buzzes are attempted captures, not successes, prey individuals or moth-specific mortality. Passes are acoustic detections, not individual bat counts. Thus these results do not quantify wetland predation pressure or a pathway to plant decline.

Keep as one supplementary experimental programme pending final eligibility and overlap adjudication, with zero additions to the main catalogue. Prioritize the per-pass endpoint when addressing activity-normalized attempts; do not combine it with count endpoints as if they measured the same ecological effect. Reproduce with `python database/alan_experiments/zeale_sources/extract_feeding.py`.

Independent calculation from raw rows reproduced all12 point estimates, variances and intervals; the covariance matrix was positive semidefinite within numerical tolerance and rank7. Independent review also confirmed ratio rounding and endpoint support.
