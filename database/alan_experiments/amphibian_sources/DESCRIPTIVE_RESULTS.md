# Shared-light-group descriptive outcomes, 2026-11-20

Source: [Dananay & Benard 2018 original workbook](https://doi.org/10.5061/dryad.pq54gh5), checksum-verified Zenodo mirror4932781; methods evidence and unresolved eligibility are in [review](../search/amphibian2018_review.md). This is a descriptive extraction, not a new meta-analysis or adopted effect catalogue entry.

## Estimand fixed for this extraction

Each source-labelled spatial group contains four mesocosms with the complete lid-by-toad treatment combination. Average these four equally, then average the four illuminated and five natural-light groups separately. This is an average across that experimental treatment mixture, not an unrestricted wetland-community effect. The supplementary design states lamps were centred within illuminated groups; workbook Metadata defines Block as a spatial grouping of four mesocosms. This supports group-level descriptive aggregation, but no mesocosm-ID physical lamp diagram was recovered. Original subgroup values remain in the archived workbook.

Use the last recorded algal date (Leaf2July, Water2July), not an average of June and July. PredatorTotal and NonPredatorTotal are totals across both pipe-sampling days according to Metadata; they are not final-time snapshots. Retain endpoint names: the NonPredatorTotal metadata confusingly says 'non-amphibian predators', so no independently validated trophic classification is claimed.

|Source endpoint|Illuminated group mean|Natural-light group mean|Descriptive ratio|
|---|---:|---:|---:|
|Leaf2July|49.1875|79.79775|0.616402|
|Water2July|21.6158125|20.9569|1.031441|
|PredatorTotal|5.6875|4.55|1.25|
|NonPredatorTotal|1.5|0.85|1.764706|

Algal endpoints are the workbook chlorophyll/relative fluorescence indicators, not dry biomass or plant reproductive performance. Invertebrate counts describe colonization/sampled presence, not attacks, consumption or predation rate. Different response dates and sampling windows preclude interpreting this table as a contemporaneous mediated pathway. Lids alter light as well as colonization. Do not infer causal mediation from opposite directions of these marginal means.

## Verification and limits

`python database/alan_experiments/amphibian_sources/summarize_groups.py` writes36 group-endpoint rows and4 arm summaries. It verifies source MD5, complete one-to-one36-ID joins, light agreement, four factorial combinations per group and no missing endpoint values. No missing records are converted to zero. Group counts are4/5, not16/20 independent lights; responses from the same groups are dependent. No interval, significance test, lnRR variance, main-meta adoption or independent-study count increment is made. Final source and indoor-phase review remain incomplete; an unsuccessful EuropePMC retrieval and an NCBI BioC HTTP200 error page did not supply primary text.

An independent reaggregation reproduced all four means and ratios, verified complete2x2 treatment combinations and no missing observations, and confirmed the trophic-label ambiguity. Workbook Metadata documents exclusion of mesocosms13–16 because of shade; they are not reconstructed as observed zeros. The original metadata gives no measurement units for the two chlorophyll columns, so none are assigned here.
