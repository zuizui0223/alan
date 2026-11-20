# Hong Kong paired spider/herbivory extraction — 20 November 2026

[Guo et al2025 primary paper](https://doi.org/10.1007/s00114-025-02031-w), [open fulltext](https://pmc.ncbi.nlm.nih.gov/articles/PMC12521300/), [original author data/code v1](https://doi.org/10.6084/m9.figshare.27273918.v1).

## Verified design

14 forest plots, seven lit–unlit pairs, four Hong Kong sites, April–June2024. Pair allocations: LammaIsland2, KadoorieCentre2, LittleHawaiiFalls2, LamTsuen1. Each plot100m², paired plots≥50m apart, pairs≥250m apart. A25W6500K cool-white solar LED at2m was installed in each plot; control lamps remained off. Reported47lux is at0.1m, NOT mean whole-plot exposure; light fell to0lux at4m. Random allocation of lit member is not explicitly reported in methods, so paired controlled experiment rather than verified randomized trial. Continuous on-time and exact installation dates by pair are not present in these data. Wet-season forest cannot be labelled wetland/riparian on that basis; coordinate-derived wetness is pending.

Spider observations are female Nephila pilipes counts, size and a length-sum biomass proxy, sampled approximately fortnightly. The417 spider records (201lit/216unlit) are NOT independent light replicates or necessarily unique individuals. Seven paired comparisons are the denominator. Leaf sampling used10leaves/species/plot/round;840leaves are subsamples. Published data contain plot-level means, not individual leaves.

## Extraction and check

Original Figshare zip49910106 was downloaded and hashed. Derived tables and executable extraction are in database/alan_experiments/hongkong_sources/. Final spider endpoint uses each pair’s latest common observation date; final herbivory uses round3. All seven pairs have both treatments on the same date, June21–27. Final n=7 per treatment, paired. Original count sums201/216 match paper. Missing values/zero-spider plots were not dropped: final counts include one lit zero. Ratio of overall means is well-defined and no pseudocount was added.

| Endpoint | Lit mean ± among-plot SD | Dark mean ± among-plot SD | RR | Paired approximate95%interval |
|---|---:|---:|---:|---:|
| Female spider count/100m² |2.1429±2.0354|2.8571±3.2878|0.750|0.447–1.258|
| Psychotria asiatica herbivory% |4.8820±2.5469|8.0121±9.0932|0.609|0.204–1.818|
| Aporosa herbivory% |2.9785±2.1339|3.3197±3.0064|0.897|0.319–2.523|

lnRR=log(meanlit/meandark). Sampling variance uses the measured across-pair covariance: [var(lit)/meanlit² + var(dark)/meandark² −2cov(lit,dark)/(meanlit×meandark)]/7. Intervals use delta-methodSE and t6; these are descriptive small-sample approximations, not the authors’ longitudinal model. Raw mean-difference variance is also supplied. Cross-endpoint covariance is calculated from paired influence values, retained in final_effect_covariance.csv; these three responses share the same experiment and cannot count as three independent studies. Four-site grouping may leave residual clustering beyond seven pairs; do not claim global precision from these intervals.

## Main interpretive limits

The paper’s spider longitudinal models found no detectable light effect. The extracted final count RR<1 does not establish a decrease: its interval crosses1. The herbivory headline concerns early wet season, especially Psychotria; final-round estimates reverse that direction and are uncertain. Under the already fixed final-time protocol, report these final endpoints while preserving time dependence in the narrative. Do not switch to early time merely for a positive effect.

Predicted spider-mediated cascade was not demonstrated: moth abundance/prey interception were not quantified, spider responses were undetectable, and plant herbivory did not correlate significantly with spider biomass. This is mechanistic counterevidence to a simple universal increase in spider predation under illumination, without proving no effect.

## Source version discrepancies requiring retention

1. Final paper names Aporosa dioica; original README and R comments call it Aporosa octandra; raw values identify only Aporosa. Preserve genus-level endpoint pending name/voucher clarification; do not perform EIV or species-trait joins as though identity were settled.
2. Figshare v1 predates final paper. Herbivory deposited R uses random=~1|Site_pair; final methods describe location plus sampling-round nesting. This final-time extraction bypasses those model differences but is NOT a reproduction of the final longitudinal model.
3. Final spider-size measurements are conditional on presence, and zero plots create missing mean size; no size effect was added to avoid changing the final eligible sample. Biomass proxy can be extracted later but must not be interpreted as measured mass.

No master tables were edited and no pooled analysis was run. Three candidate final-time effects now have raw-unit provenance and paired covariance; eligibility and taxonomy flags accompany them.
