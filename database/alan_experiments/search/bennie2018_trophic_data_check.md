# Bennie2018 trophic data: independent check

## Primary design and endpoint

[Methods](https://doi.org/10.1111/1365-2664.13240) describe54 Penryn mesocosms, six randomized units per lighting×trophic cell. Thirty-six received herbivores. Slugs were counted monthly for five minutes under standardized survey lighting; experimental lights were switched off on survey nights. Predator abundance was not systematically measured. Slug analysis stops at September2014 to avoid later reintroduction. Published models use monthly counts and mesocosm random effects; an equal-month per-plot average is a different estimand. Plot letters are not reported as paired blocks. Counts do not measure individual attack or consumption rates. Source: cached `bennie2018_trophic/primary.txt`, Methods and Results; [original archive](https://doi.org/10.5061/dryad.8217n76).

## Independent numeric checks

Inspected original `slugs2013.csv` and `slugs2014.csv`, separately from the extraction script. Both contain the same36 unique plots and treatment labels; each of six cells contains6 plots. Six2013 columns (June–November) plus nine2014 columns (January–September) give540 observed plot-month values, without missing values. December2013 is not supplied. The summary averages15 observed months, not a continuous16-month exposure-weighted series.

| Trophic group | Light | Per-plot period mean, then arm mean | September2014 mean |
|---|---|---:|---:|
|Grazer|Control|6.900000|0|
|Grazer|Amber|6.266667|0.333333|
|Grazer|White|6.244444|0.333333|
|Predator|Control|4.911111|0|
|Predator|Amber|4.955556|0.666667|
|Predator|White|2.988889|0|

Both final dark groups have all six observations zero. All four final light/dark lnRRs are undefined; no pseudocount, infinite finite effect, or substitute final month is justified. Independent n remains6 per cell, not90 monthly observations.

For the supplementary period average, scalar calculations using arm means and sample variance/6 reproduced all four lnRRs and sampling variances in `bennie_trophic_sources/slug_period_effects.csv`. Predator-minus-grazer contrasts reproduced:

| Light | lnRR difference | Variance | Ratio of response ratios | Approximate t5 interval |
|---|---:|---:|---:|---|
|Amber|0.1052859004|0.0909555460|1.1110282083|0.511723–2.412210|
|White|-0.3967692704|0.1707577148|0.6724891709|0.232467–1.945401|

The shared-control covariance between these contrasts is0.0543750082, independently obtained by adding the two dark-cell relative mean variances. Interval endpoints match the saved t5 calculation. These are secondary descriptive approximations, not a reproduction of the original longitudinal GLMM.

## Same-programme check and historical discrepancy

The recovered `PeaAphid2013.csv` joins existing `bennie_sources/repeated_original_units.csv` on all378 plot-date keys with matching treatment identities. Final12September2013 values match for all54 plots (including corresponding missing entries). Across the complete time series, two earlier values differ:

| Plot | Date | Recovered2018 archive | Existing2015 source extraction |
|---|---|---:|---:|
|B6|2013-08-29|7|0|
|C4|2013-08-29|4|0|

Do not describe the whole2013 time series as numerically identical or silently overwrite either source. This discrepancy does not change the currently retained final-date effect. It does not establish a new independent programme. No source data were modified in this check.
