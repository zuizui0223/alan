# Toad feeding experiment: recovered contrast and calendar-confounding audit

Rechecked2026-11-21. [Primary article](https://doi.org/10.1038/s41598-020-63503-9); Supplementary TableS2 pp4–5. The existing source extraction was recovered from the repository and checked against the rendered supplement; this is not a newly acquired study. A fresh official Europe PMC supplement download matched the stored PDF hash.

## Direct feeding endpoint and descriptive contrast

Six outdoor4×4m enclosures near Townsville were used repeatedly, with five newly introduced toads planned each night. The Results report171 dissected toads across37 enclosure-night trials; most sites have6 trials, site4 has7. The response in TableS2 is mean gut-content mass per toad per overnight trial. Individuals are subsamples, not171 independent lighting units. Lamps were approximately7.3–7.8lux when on; ambient off-state values ranged0–3.65lux, so the control was not uniformly dark.

The existing calculation averages nights within each enclosure×treatment, then weights the six enclosures equally. Mean gut content was2.549g under lighting and0.154g when off, ratio16.574. One enclosure's off-state mean is zero, but the across-enclosure mean remains positive; no pseudo-count is used. The average ratio is a ratio of means, not a mean of six ratios. This is overnight food mass retained in guts, not an observed attack rate, prey population loss, or quantified adult-moth predation. The paper reports a greater fraction of flying taxa in lit gut contents, but the extracted TableS2 does not identify prey taxa per trial.

## Why the existing interval is not a pooled causal-effect estimate

`calendar_treatment_support.csv` shows14 distinct month/day combinations. **Every observed date contains only on or only off trials.** Sites4–6 were sampled in January and sites1–3 in March/April. Thus all sampled enclosures on a date share the treatment, and there is no simultaneous within-date lit/dark contrast. Adding unrestricted date fixed effects would make the light indicator redundant; it does not identify an adjusted light coefficient.

The previously saved t5 delta-method interval5.435–50.541 assumes six independent enclosure-pair contributions. It retains within-enclosure covariance but does not account for shocks shared by enclosures sampled on the same date or remove temporal confounding. It is explicitly labelled a **naive descriptive interval, not adopted for pooling** in `effect_paired.csv`. Neither this interval nor the37 trials is used to increase the precision of an independent predation meta-analysis. A large point ratio is not a solution to the dependence problem.

Original author models adjust measured environmental covariates and use an enclosure random effect; the extracted marginal contrast is a different estimand and is not a replay of that adjusted model. Ambient light and lunar phase interactions in the original paper should not be inferred from this single marginal ratio. The available calendar does not prove random assignment of on/off nights, and the extra off night at site4 is retained rather than removed to enforce balance.

## Source handling and manuscript use

`trials_primary.csv` preserves all37 reported trial rows; `enclosure_treatment_means.csv` preserves all12 cell means. `extract_primary.py` reproduces the original source extraction, point ratio and historical naive interval and now records the calendar limitation. TableS2 gives day/month; year2018 is inferred from the original December2017–April2018 study window, not separately printed in that table.

The enclosure grassland/forest setting is not a wetland treatment merely because the predator is an amphibian. This experiment supplies direct feeding evidence with a temporal-confounding limitation, not proof of the pathway from water reflection through predator attacks to wet-plant decline. It is presented separately in supplementary Results, with no change to the108-row effect catalogue or114-row reported-statistics registry.
