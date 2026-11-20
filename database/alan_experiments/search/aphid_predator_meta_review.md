# Aphid predator-modification synthesis: independent audit

## Primary endpoint/design compatibility

[Penryn first report](https://doi.org/10.1098/rstb.2014.0131) measures Acyrthosiphon pisum with standardized timed searches; six randomized mesocosms per lighting×trophic cell contribute. Its retained12September2013 endpoint is final for that report, not terminal for the programme, which has later follow-up. Predator treatments combine ground beetles and ladybirds. [Miller2017](https://doi.org/10.1098/rspb.2017.1195) counts the same aphid species after9days using whole-cage collection; four blocks contribute to each predator/no-predator lighting comparison. The complete unheated subset has24 cages; warmed allocation errors remain excluded and unresolved. Adjacent-row blocks are reconstructed from the stated layout.

These are distinct experiments with different host vegetation, duration and predator menus. Both support experimental effect modification of recorded aphid counts, not direct consumption. Timed-search versus census ratios are comparable only under treatment-relative detection assumptions; logarithms do not remove differential detectability. Miller Methods report solar-powered LED lamps but do not establish white colour or spectral composition.

## Independent numerical verification

Recomputed Penryn from original `PeaAphid2013.csv`, using the final date and raw treatment cells. All36 relevant records have the expected six plots per cell. Computed light×predator differences directly from four cell means; variances are sums of independent cell relative mean variances. The off-diagonal covariance is the sum contributed by the two shared dark cells. No plot-letter pairing was introduced.

Recomputed Miller from the original `Miller data ProcB.xlsx`, sheet `Field data`, restricting to Warm=N. Reconstructed four blocks, then formed block-level influence scalars for each predator-minus-none contrast. Their sample covariance divided by4 reproduces the saved2×2 covariance, including cross-treatment block dependence and the common no-predator comparator.

Both two-contrast inputs and their full covariance matrices agree. Within-programme equal averaging gives:

| Programme | Mean log modification | Sampling variance |
|---|---:|---:|
|Penryn first-report endpoint|−0.10364468934|1.16186997956|
|Miller unheated final endpoint|−0.50722865006|0.25497811769|

The first estimate averages two light colours; the second averages two predator identities. They are fixed experimental-menu averages, not representative random samples of colours/predators. Within-programme variance uses all covariance entries divided by4; dropping the shared terms would materially understate uncertainty.

Independently minimized the restricted likelihood for the sole between-programme residual difference, rather than reusing the implementation's general weighted-mean likelihood. REML heterogeneity is at the zero boundary. Recomputed modified Hartung–Knapp max(1,q), df1:

- pooled log modification:−0.434599072593;
- response-ratio ratio:0.647524225251;
- 95% interval:0.001940693302–216.050429901;
- q=0.114959404374; modified variance=0.209091871572;
- weights:Penryn0.17996151, Miller0.82003849.

Saved [pooled output](../aphid_predator_synthesis/pooled_estimate.csv) agrees to relative tolerance1e-10. Numerical verification passes for [synthesize_aphid_predator.py](../synthesize_aphid_predator.py).

## Interpretation and remaining qualification

A ratio below1 means lighting's proportional effect on recorded aphids is smaller with predators under the specified menus. It does not measure the fraction eaten, prove a mediating pathway, or estimate plant reproductive consequences. The interval includes1 very widely; neither a general amplifying interaction nor absence of interaction is established. Zero estimated heterogeneity with two programmes does not imply uniform biological responses.

The complete unheated subset and first-report endpoint rules must be stated before reuse; do not describe either as all treatments or terminal programme data. Keep the original author longitudinal/factorial findings distinct from these report-specific endpoint calculations. The slug-period, plant-biomass and other invertebrate endpoints are not part of this synthesis. Common-sample dependence prevents counting Penryn's later publications as independent corroboration.

Source files: [Penryn2018 recovered archive](https://doi.org/10.5061/dryad.8217n76) (its final2013 aphid values were previously matched to the first-report source), [Miller original archive](https://doi.org/10.5061/dryad.h7k26). The earlier two Penryn source-version discrepancies are on29August2013 and do not affect the selected12September endpoint. This audit changes no source values or eligibility flags.

Implementation follow-up: the unsupported Miller `white_LED` label was replaced
with `LED_colour_unreported`, and the Penryn endpoint was explicitly labelled
as final in the first report, not terminal in the programme. Numerical values
were unchanged by these metadata corrections.
