# Zhou et al. 2023: manipulated gradient, not replicated plant light/dark contrast

Read2026-11-20: [preprint](https://doi.org/10.32942/x2b60s) and [published article](https://doi.org/10.3390/ecologies4040046), both Methods/Results. PDF source hashes are in source_manifest.json. Independent reading confirms the design and version assessment below.

## One experimental programme, two reports

Same ordered authors, species, XTBG location (21d54mN,101d16mE), November2019 installation/November2021 leaf sampling, LMA slopes/intervals and soil nutrient ratios identify the same experiment. Prefer final version; do not count the preprint as a separate experiment. Final version adds sample sizes and conceptual distant-control explanation. Final DOI is a supplemental report outside the current indexed search, so the search-only publication-family table is not used to pretend that it was an indexed hit.

## Measurements and eligible comparison

All five plots received a10W LED; leaf traits were sampled in two plots selected for natural plant distribution, one plot per species. Final Methods report35 Colocasia gigantea individuals/165leaves and19 Melastoma candidum individuals/95leaves. Individual plants and leaf discs are not independent illumination replicates. Lights operated20:00–05:00, at2m height; final source reports6696K. Plant models use transformed lamp distance, canopy openness and interaction, with individual random intercepts. Distant plants are described as conceptual controls, not independently assigned unlit plots. Keep as a manipulated spatial gradient; do not bin distance after seeing outcomes to create a two-arm lnRR. One species per plot confounds species or shade-tolerance differences with plot differences.

Soil C/N/P were sampled before installation (June2019) and afterwards (June2022), at lamp and10m positions in five plots, three subsamples per location. This is a separate paired spatial-before/after estimand with five potential plot units, not35/19 plant replicates. The June2022 soil measurements do not precede the November2021 leaf endpoint; do not claim measured temporal mediation.

## Results and mechanism boundary

Final Table1 reports ALAN-distance predictor coefficients for LMA: C.gigantea -0.1052,95%credible interval[-0.1500,-0.0613]; M.candidum -0.0434,[-0.1147,0.0278]. These are model coefficients for transformed predictors, not log response ratios or effects per lux. Their units/scaling must be inspected in original code before quantitative synthesis.

Soil nutrients declined in both near and distant samples; in the C.gigantea plot some declines were smaller near lamps. Thus 'nutrient increase' would misstate the absolute time trend. The published Methods/Results do not report measured insect carcass flux, nutrient transfer from insects, pollination, predation or plant population decline. The archive caveat below prevents concluding that no insect data were ever collected. Insect-mediated fertilization is the paper's proposed mechanism, not an experimentally isolated effect. A cited ant study used the same experimental setup; any later ant/plant synthesis must check shared plots rather than count it as independent evidence.

Decision: hold outside the binary light/dark plant meta; retain manipulated-gradient trait and soil context. No added main-catalogue effect and no claim that two species establish contrasting functional-group responses.

## Archived original data/code recovered

Zenodo concept record10060043 resolves to version10060044; archive `congz007/light_project-light_project_v0.1.zip` (11,800,312bytes) was downloaded and its MD5 verified against deposited metadata. It contains repository snapshotcc45e9f. Plant rows reproduce165 leaves/35individuals/1plot for Colocasia and95/19/1 for Melastoma. `soil_alan.csv` has45 rows across distance positions, with missing nutrient values at intermediate distances; these must not be treated as45 replicate nutrient observations.

`R/stan_clean.R` reveals the ALAN predictor is within-species scaled log(1/distance^2), with scaled canopy openness and their interaction; response is log(LMA). Thus Table1 coefficients are changes in log(LMA) per one sample-standard-deviation of that transformed distance proxy at centred canopy openness, not per lux or a light/dark response ratio. Stan includes individual random intercepts; it does not supply replicated species-by-plot contrasts. MCMC was not rerun and posterior outputs were not independently reproduced.

A crucial caveat: `scripts/Insects.R` references `data/Drymass_funnel.csv` and sums two time labels before plotting drymass by distance, but that CSV is absent from this frozen archive. The script alone is not evidence of measured carcass deposition or nutrient transfer, yet means that 'no insect data were collected' would be unjustified. A publicly available raw insect deposit or explanation of its sampling method is a new acquisition lead. Do not execute its dry<0.5 plotting filter as a meta-analysis exclusion rule. No insect effect or mediation coefficient can currently be extracted.

## Follow-up: current repository and related ant report

On2026-11-20 the public GitHub default-branch tree was complete (not truncated), at commitcc45e9f18bb4f683bd3ca24ef038204e26325837, identical to the archived release snapshot. `data/Drymass_funnel.csv` was absent. The GitHub commits endpoint restricted to that exact path returned an empty list. This rules out recovery from the checked current tree/exact-path default-branch history, not from all possible renamed files, other branches or external deposits. Stop repeating that same retrieval without a new lead.

Reference11 is Yang & Nakamura2022, [Cavity-dwelling ants tend to colonize close to artificial light](https://doi.org/10.17520/biods.2022067). Publisher abstract is available via [indexed publisher record](https://www.biodiversity-science.net/EN/abstract/article/1005-0094/82308); direct article and linkedPDF requests returned403. Do not count it as full-text read. Abstract describes15 lamps across primary forest, secondary forest and rubber plantation; bamboo nests at0,10,50,100m for seven weeks in dry/rainy seasons. Increased occupancy beside lamps was reported only in dry season, not a universal seasonal response. Ant community composition did not change with distance according to the abstract.

This offers a predator-location/colonization endpoint, not measured predation pressure or insect-derived nutrient flux. The plant paper says the same setup was used, but plot identifiers, deployment dates and exact overlap still require full methods/data. Do not combine the plant and ant reports as independent programmes, or claim simultaneous plant-ant-soil mediation. Record as a citation-traced lead outside the indexed search; no effect extraction or search-hit inflation.
