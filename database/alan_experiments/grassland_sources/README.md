# Lizard grassland: final2015 dry-material extraction

[Primary Methods/Results review](../search/grassland2018_review.md); [original archive](https://api.figshare.com/v2/articles/29731352). Archive download, hashes and all15 members are recorded alongside this file.

## Source and units

The biomass workbook has120 plot-years,30 plots and five invariant treatments across2012–2015. Each final-year arm contains six plots. The paper describes random allocation, not matched blocks; plot-name letters are not used to invent pairing. All2015 totals equal their component sums. The single earlier mismatch (G7,2012) remains in a separate table and does not enter the final-year estimates.

The source `total` includes `Dead`. We retain this endpoint and the dead-material endpoint, and derive `live_total` by summing the remaining named plant/material columns. These are overlapping quantities, not independent confirmation. Viability was not measured beyond the author's material separation.

An undated harvest CSV scales matched2012 values by four. Absolute area-normalized units are therefore not assigned here. Any common multiplicative conversion within an endpoint cancels from lnRR and its covariance; raw means retain `source_workbook_mass_units`. Final2015 has no duplicate harvest table to append.

## Quantitative results

All68 contrasts (17 endpoints×four light arms) are retained.28 have positive means in both arms;40 have a zero mean and no finite lnRR. No pseudocounts, abundance threshold or significance filter are used. Covariance accounts for four comparisons sharing the same dark plots and multiple measurements on each plot; its rank is25. A meta-analysis must select a defined endpoint, not fit all overlapping totals and parts as independent observations.

For amber/high-white/low-white/part-night-white respectively, live-component response ratios are0.616,0.709,0.719,0.746; totals including dead material are0.706,0.728,0.747,0.763. The latter four approximate Welch-delta intervals all include one. Amber live-component ratio has an unadjusted interval0.417–0.910; this is one exploratory contrast among correlated comparisons, not a multiplicity-adjusted general conclusion.

Agrostis ratios are below one in all arms, while Anthoxanthum and Holcus ratios are above one; several intervals are very wide. These are final-year estimates, not reproduction of the original repeated-year models or a measured pollination/predation mechanism. Rare-species zero-arm comparisons remain visible and cannot be silently dropped from claims about the community.

Run `python database/alan_experiments/grassland_sources/extract_final_biomass.py`. Finite entries enter the provisional catalogue as extracted, not finally adjudicated.
