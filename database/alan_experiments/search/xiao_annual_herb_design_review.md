# Xiao annual-herb field experiment: design review

Reviewed 2026-11-20. Sources: [preprint DOI](https://doi.org/10.1101/2022.12.11.519667), [v2 source XML](https://www.biorxiv.org/content/early/2025/08/21/2022.12.11.519667.source.xml), [final article DOI](https://doi.org/10.1016/j.biocon.2025.111472). Independently read Methods and Results from the recovered v2 XML, and accessed v1 PDF via requests (HTTP 200; 947398 bytes). Publisher full/PDF access through the web tool returned 403; final-paper Methods, tables and supplements were not independently verified. The final abstract/title were indexed and accessible. This note does not claim a final-full-text review.

## Source-derived summary (under 200 words)

The 2021 Inner Mongolia experiment used 20 open-cage plots within a 40-plot light-by-animal-exclusion programme. Ten plots received 6500-K ALAN (21.15 lux); ten had dummy fixtures. Ten marked Elsholtzia densa plants per plot were monitored daily until September 23. Biomass sampling used additional individuals in August. Plot identity entered individual-level models.

V2 reports earlier budding, flowering, fruiting and seed maturity, with flowering advanced 6.2 days. Fruiting inflorescences decreased 25.2%, but support met only the authors' 90%, not 95%, credible-interval criterion. Inflorescence biomass fraction and total biomass showed no significant treatment response. Fruiting-inflorescence count was explicitly a reproductive-capacity proxy: seed output was not directly quantified because seeds shed prematurely. Pollination disruption was proposed, not measured. V2 states that data are available on request.

## Own eligibility and extraction audit

- Lighting replication is **10 plots per arm**, not 100 marked individuals per arm or 45 versus 52 harvested individuals. The three-to-six biomass plants per plot are unequal subsamples. Estimate plot-level means or explicitly model the nesting; retain common-plot covariance across traits and dates.
- Only open-cage plots contribute here. This is an ALAN comparison under animal access, not an ALAN-by-pollinator-exclusion test. No measured pollen deposition, single-visit efficiency, pollen supplementation, visitation frequency, predator attack or plant-animal phenological mismatch can be inferred.
- A final reproductive endpoint candidate is **number of seed-bearing inflorescences per marked plant**, with plot aggregation. It is not fruit-set probability, seeds per flower, viable-seed production, recruitment, or lifetime fitness. Inflorescence length and undeveloped-top-inflorescence frequency are separate morphological/developmental outcomes.
- Phenology is naturally expressed as a difference in days; a ratio of Julian dates is origin-dependent and not suitable as a biological lnRR. Repeated phenophases share individuals/plots. The August biomass harvest and September reproductive census represent different times and subsamples.
- Random selection of plants or redistribution of seeds is not evidence of random lighting assignment. The inspected Methods do not explicitly establish the allocation randomization. Retain this uncertainty until a layout/protocol resolves it.
- The latest preprint uses Bayesian model-specific likelihoods and transformed responses. Its estimated coefficients/credible intervals cannot automatically replace arm means, SDs and independent-plot sample sizes. The marginal fruiting result must not be strengthened by copying the broader abstract wording or older-version claims.

## Provenance and decision

The preprint and final publication are one study lineage, not two independent experiments. Keep version differences explicit and prefer verified final values when available. This field programme is separate from [Xiao2025 greenhouse herbivory](https://doi.org/10.1098/rspb.2024.2811), which used multiple plant species and Spodoptera larvae in Changchun. Do not import its compartment IDs, species rows or Dryad files here.

The field study belongs to a larger long-term grassland platform; any later reports from that platform require campaign/plot/control overlap checks before counting independent programmes. No public raw-data repository was identified in the inspected v2 data statement, and no raw observations were extracted. Request-available is not evidence that data do not exist. Verdict: eligible field reproduction-proxy/phenology candidate; quantitative extraction pending plot-level source data and verification against the final paper. It cannot independently demonstrate a pollination mechanism.

## Independent Table 2 transcription check

Visually inspected rendered v2 PDF page 32 (`table2.png`) on 2026-11-20. All five requested ALAN coefficients and intervals agree with the table:

| Response | Estimate | 95% credible interval | 90% credible interval |
|---|---:|---|---|
| Budding date | −3.588 | −4.954 to −2.243 | −4.713 to −2.464 |
| Blooming date | −5.494 | −7.673 to −3.457 | −7.294 to −3.785 |
| Fruiting date | −3.335 | −5.107 to −1.583 | −4.814 to −1.885 |
| Seed maturity date | −4.526 | −6.609 to −2.462 | −6.259 to −2.784 |
| Fruity inflorescence number | −0.285 | −0.628 to 0.050 | −0.567 to −0.008 |

These are author-model estimates, not independently calculated arm effects. The phenophase models are Gaussian, while inflorescence number uses a negative-binomial model. Its 95% interval includes zero and 90% interval excludes zero. The text's 6.2-day flowering advance differs from the model's 5.494-day contrast; retain both with their source/estimand labels rather than substituting one for the other. Table sample sizes of 200 refer to individuals, not independent lighting plots. This check did not rerun models, recover raw means/SDs, or infer SEs from credible intervals.
