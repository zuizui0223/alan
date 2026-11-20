# Supplemental pollination search: remove mandatory experiment terms

Executed2026-11-20 through EuropePMC, with date restriction through that date. Exact query, raw-response hash and counts are in `pollination_expansion_manifest.json`; all51 returned records are in `pollination_expansion_records.csv`. The raw response and abstracts are cached externally. This additional search omits the mandatory experiment/manipulation concept, because field studies can describe their designs without those words. It supplements—not replaces—the521-record search, citation tracing, and unfinished WoS/Scopus searches. It is not a complete new systematic review.

Exact-ID/DOI reconciliation finds26 original records and one already supplementary-registered report (Knop2017). The remaining24 records require identity and title/abstract review; they are not24 new eligible studies. Missing DOIs, preprints and duplicate index sources require further reconciliation. No total PRISMA count or eligibility decision was silently overwritten.

## Priority leads from abstracts, not extracted effects

- **10.1093/icb/icab010**, Direct and Ambient Light Pollution Alters Recruitment for a Diurnal Plant-Pollinator System. [Primary index](https://europepmc.org/article/MED/33724371). Yucca–Tegeticula mutualism; abstract describes field light towers with dark trials for adult captures, and separate fruit/larval recruitment measurements. Direct-light manipulation and observational skyglow must be separated. Retrieve full Methods and original data before identifying independent units or pooling. Positive recruitment associations must not be excluded for conflicting with a decline hypothesis. Not added to primary-read counts.
- **10.1007/s00114-021-01783-5**, Negative effects of light pollution on pollinator visits are outweighed by positive effects on the reproductive success of a bat-pollinated tree. [Primary index](https://europepmc.org/article/MED/34994866). Abstract describes contrasting urban light/noise sites and structural-equation modelling of Ceiba pentandra, not randomized illumination. Retain as an observational mechanism lead, separate from experimental pooling; reported visit and reproductive responses differ in direction. No causal effect or numerical estimate extracted, full-text classification pending.

Neither abstract establishes a wet-versus-dry effect or a mechanism for the macro wet-affinity result. The remaining records are not yet screened in this new search; historical decisions are referenced through exact matches rather than copied as final eligibility.

## Boom2020 supplement retrieval

The previously registered Boom2020 supplement remained unrecovered: a second EuropePMC supplementaryFiles request timed out after55seconds. A publisher-CDN attempt failed DNS resolution. No supplement contents, covariance or raw-data availability conclusion was inferred from these failures. The primary Table1 extraction remains valid, but its original variance holds remain.
