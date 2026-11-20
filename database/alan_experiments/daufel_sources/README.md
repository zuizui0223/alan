# Conditional oviposition effects, Daufel2023

Primary [Methods/Results](https://doi.org/10.1111/eea.13346) describe eight
paired outdoor blocks, each with one shielded3000K lamp and lit/dark pools.
Measured illumination was13.46 versus0.24lux; lights operated throughout the
experiment. This is not a night-only on/off exposure. Pools were2m apart and
blocks3m apart, allowing redistribution of oviposition among nearby pools.
This measures deposition/site selection, not total population reproduction,
emergence, predation or pollination.

Species assignment is partial:228/261 rafts were reared and identified as
Culex restuans; only11 individual eggs were identified as Ochlerotatus japonicus.
Keep endpoints named by deposition form rather than claiming every egg's
species was verified. The paper reports261 rafts and197 individual eggs over
seven days. Raw-data discrepancies below remain unresolved.

## Original-source extraction

[Figshare19230045](https://doi.org/10.6084/m9.figshare.19230045) supplies the
workbook; published MD5 is verified. `audit.json` records its URL and hashes.
Numeric Pool AND Pair identifiers select16 data rows per sheet, avoiding
numeric summary rows. `extract.py` creates paired-unit data, missing-day cells,
two conditional effects and their2×2 sampling covariance.

All16 pools have six daily entries within seven dated columns. Pair1,2,7,8
lack7June; pairs3–6 lack3June, identically for both arms and endpoints.
Stored cumulative values exactly equal the six recorded entries. Missing
cells are not filled or interpreted as zero; no complete seven-day rate is
claimed. Counts sum to261 rafts and193 single eggs, the latter differing from
the paper by four. Species-labelled workbook sheets are not substituted for
the two deposition sheets.

| Recorded cumulative endpoint | Lit mean | Dark mean | Ratio | Approximate95% interval |
|---|---:|---:|---:|---:|
| Egg rafts |10.625|22.000|0.482955|0.273270–0.853535|
| Single eggs |12.500|11.625|1.075269|0.303673–3.807392|

Delta-method variance uses the eight paired influence scores
`lit/mean(lit) − dark/mean(dark)`; interval quantile is t7. Zeros remain in
arithmetic means. This is not a re-fit of the original Poisson mixed model.
All eight raft contrasts are negative. A second analyst agent independently
recomputed both estimates and covariance; this is not a two-human extraction.

No new catalogue or pooled effect is adopted pending temporal completeness,
single-egg discrepancy and programme adjudication. The conditional raft
result can be shown separately in a provisional Results section. Do not add
these endpoints to adult light-trap capture or predator meta-analyses.
