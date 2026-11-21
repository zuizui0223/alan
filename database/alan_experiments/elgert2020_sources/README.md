# Glow-worm mate attraction: field test separated from laboratory behaviour

2026-11-21. [Elgert et al. 2020 primary article](https://doi.org/10.1098/rspb.2020.0806); [author data and Methods](https://doi.org/10.5061/dryad.2z34tmphw). This extends the existing batch 2 review and does not add a new paper to screening counts.

## Field design

The study used six shared sites around Tvärminne, Finland, from 6 June to 7 July 2017. A pole-mounted white LED was on for 38 trials and off for 19 trials. Every trial had two dummy female traps of equal glow intensity: inner position at 1.1 m from the pole, 15–20 lux when lit; outer position at 2.3 m, 1–1.5 lux when lit. Three lure intensities varied among trials. Both positions share the same lighting allocation and night. The 114 trap observations are not 114 independent lighting replicates, and the six sites are not six lit plus six separate control sites.

The outcome was at least one male caught per dummy trap, not male counts or completed mating. The author binomial-logit GLMM included light, position, their interaction, and initially lure intensity, with site nested within date as random terms. Lure intensity was removed from the final model. The article states equal treatment distribution across sites, but the detailed night-by-site allocation cannot be reconstructed from the totals alone.

## Reported quantitative evidence

`reported_field_tests.csv` preserves all three Table 2 tests: light F(1,110)=6.489, P=0.012; position F(1,110)=7.683, P=0.007; light × position F(1,110)=5.767, P=0.018. The narrative states that male attraction under lighting was lower in the inner than outer position. This is a position-dependent lighting response, not a common lighting effect independent of distance. Position changes distance and illumination together, so this design does not isolate their causal contributions or establish a universal lux threshold.

The table supplies tests, not signed coefficients, group probabilities, SEs or coefficient covariance. No odds ratio or sampling variance is inferred from F or P. The numerator/denominator degrees of freedom are reproduced as printed, not converted to independent site n. The outer position still receives artificial light; it is not substituted for an unlit control.

## Separate indoor experiment and environmental limits

The same paper's 26 lit and 63 control females were tested individually in laboratory arenas. The laboratory glowing counts (6/26 versus 61/63) are not the field mate-attraction outcome and are not pooled with field captures. The endpoint is neither plant pollination nor recruitment/population decline.

Natural night-light measurements of 0.1–0.6 lux were taken in 2020, whereas the field trial occurred in 2017. They must not be treated as simultaneous trial baseline measurements. The regional coordinates describe the collection/station area, not individual field poles; no experimental water distance, site wetness or island/mainland effect is inferred.

## Access and current use

The cached primary XML Methods/Results and the public Dryad Methods were read. The Dryad landing page links `Elgert_data_ESM.xlsx` (laboratory and field sheets); clicking its ordinary public download returned HTTP 403. No access restriction was bypassed. The original workbook, paired capture outcomes, coefficients and covariance remain unavailable in this extraction. The XML hash and attempted download are in `audit.json`.

This source contributes a reported field interaction to the supplementary Results, but no pooled effect. Its three tests remain in this source table; the central 114-row reported registry and 108-row group-mean effect catalogue are unchanged. The quantitative pooling gate requires trial-level outcomes and site/date allocation, or adequate signed model contrasts and covariance.
