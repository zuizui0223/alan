# Komine et al2020 toad food-intake extraction

[Primary paper](https://doi.org/10.1038/s41598-020-63503-9); original SupplementaryTableS2 downloaded through [EuropePMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7162902/supplementaryFiles). No external meta-analysis data used. PDF text extraction, not figure digitization; all37rows and12columns inspected against rendered original pages4–5.

## Design and exact endpoint

Six4×4m outdoor enclosures around Townsville,5–15km apart, in open eucalypt forest with grassy understorey. A fluorescent lamp1.5m high provided7.30–7.81lux when on. Ambient off levels differed by site:0,0.01,0.22,0,0.01,3.65lux. Therefore the sixth control is low light, not dark. Toads were fasted40hours, five randomly selected animals put into each enclosure overnight, then collected and dissected next morning. Different toads each night. Endpoint is mean gut prey mass in grams PER TOAD in that enclosure-night, not total prey removal/16m², predator abundance, attack probability, or moth-specific mortality. Variable recapture totals171toads are given only by site; no per-trial animal denominator is supplied. Do not multiply trial means by five to invent total removals.

Treatment alternated across nights.37trials: six per enclosure except site4seven. Original table actually includes19off/18on trials; site4has4off/3on, all others3/3. Sites4–6 were sampled in January, sites1–3 in March/April. Year2018 inferred from main methods December2017–April2018; table itself gives day/month only. Within each cohort light conditions occur on common dates, so on/off comparisons are not simultaneous and are potentially confounded by date/weather/moon. Allocation to light sequence was not reported randomized. Site4starts with two successive off dates. This is an outdoor crossover feeding experiment, not6independent lit plus6different darkenclosures.

## Preserved observation window

Every trial is a new one-night exposure with different animals; there is no repeatedly exposed population with a single final development date. Average all reported trial-night means within enclosure and treatment, then compare the SAME six enclosures. This defines an exposure-average estimate before seeing significance and retains all37trials, including zeros. There is no exact same-date matchedlit/dark trial within each enclosure; pairing is by enclosure acrossits trialwindow. Selecting only lastnight per treatment would discard repeated exposure replicates and substitute different calendar windows. Do not label this a final-time longitudinal estimate.

## Result and uncertainty

Across six enclosure-treatment averages: lit2.549389g/toad/night (SD2.330005); off0.153819g/toad/night (SD0.162745). Measured paired covariance0.161948g². lnRR2.807829, paired delta-method variance0.188125, RR16.5739. Approximate t5interval5.435–50.541. Mean paired difference2.395569g, variance0.855252g². Site3has meanoff0 but the across-enclosure mean ispositive; no pseudocount or deletion was used. Smalln, lowbaseline and heterogeneous effects limit ratio precision.

These are DESCRIPTIVE paired enclosure contrasts, not a reproduction of the authors' weather/moon-adjusted mixedmodel and not a fully temporally independent causal estimate. The covariance formula is [s²lit/meanlit²+s²off/meanoff²−2cov(lit,off)/(meanlit×meanoff)]/6. Shared samplingdates among threeenclosures can induce residual dependence; confidenceintervals cannot account for unmeasured temporal confounding. Treat as separate crossover design stratum if used; do not silently pool as parallel randomized lit/dark trial.

## Ecological interpretation

This study directly measures increased predator food intake under artificial illumination. It does not measure insect influx independently, proportional prey mortality, long-term predator population change, pollination or plant demography. Ambient light and moon context modify feeding in the original analysis. The experiment is terrestrial openforest, not automatically wetland because the predator is amphibian. Coordinate-derived water distance and waterfrequency remain uncomputed; species ecology is not experimental habitat classification.

## Files and verification

`database/alan_experiments/toad_sources/trials_primary.csv`:37 original table rows including weather/moon. `enclosure_treatment_means.csv`:12rows (6pairedunits). `effect_paired.csv`:one descriptiveeffect. `extract_primary.py`:reproduces table parsing, uniqueness,37trials,replicatecounts and effectcalculation. `source_manifest.json`:original PDFhash/sourceURL. SourcePDF and renderedpagechecks are in externalcache toad2020. Script run successfully; no mainDBtable changed and no pooledmeta-analysis run.
