# Pollination extraction reconciliation, 2026-11-20

The screening sheet previously retained stale pre-extraction eligibility labels for Giavi and Lu. Their source-derived effect files and the catalogue already contain paired daytime effects. Updated `pollination_search.csv` records current availability explicitly; its historical `new_effects_extracted` column refers to the original screening pass, not current availability. Use the new availability columns for current inventory. No new numerical effect is claimed by this reconciliation.

| Report | Raw catalogue rows | Separate published model coefficients | Current usable comparison |
|---|---:|---:|---|
| Giavi2021 |4|0|Daytime visits to author-selected21 plant labels; three orders and their sum overlap|
| Lu2026 ECE |7|0|Daytime source-row means; total, five guilds and three-order sum overlap|
| Knop2017 |0|2|Lighting-by-all-visitor-access interactions in fruit count/proportion; two dependent outcomes|
| Macgregor2019 |0|0|Reported positive full-night effect on seed-setting success; numeric model coefficient/SE unavailable|
| Lu2026 JPE |0|0|Abstract-level reproductive findings; full Methods/variance not verified|

These are neither13 independent effects nor five independent programmes. Four plus seven raw rows agree with the current catalogue. Knop's two coefficient estimates are separate from it. Shared endpoints, controls, sites and report overlap remain explicitly relevant.

## Comparability

Giavi and Lu both allow a three-order daytime visit contrast, but selection and observation populations differ; preserve study-specific effects rather than silently broadening their target to all daytime pollination. Lu night records contain unresolved session-key and guild-sum problems. Neither daytime source measures pollen delivery or reproductive compensation.

Knop estimates a light-by-access interaction, whereas Macgregor reports a lighting main effect adjusted for access regimes. The latter must not be pooled with the former simply because both involve reproduction. Seed-setting success differs from seed numbers conditional on successful undamaged capsules. The Macgregor positive finding must remain visible even while its variance is unavailable, preventing a negative-only reproductive evidence narrative.

Macgregor source script `FinalScripts.R`, lines575–590, specifies binomial `glmer` models for full-night versus darkness with a random intercept for `fRound`; the final formula includes LitUnlit, Pollinators, Distance and usage. The binomial default is a logit link. Formula availability does not supply coefficient magnitude or variance, prove sufficient handling of rotating plots/repeated plants, or constitute a model rerun. Primary PDF Methods/Results, Figures2–4, and archived author scripts were checked. Supplemental institutional records at White Rose149121 and NERC522598 list the main article, not a recovered coefficient table.

Source links: [Giavi original data](https://doi.org/10.5281/zenodo.4540407), [Lu original supplement](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles), [Knop extraction](knop_sources/REPRODUCTIVE_INTERACTION.md), [Macgregor primary article](https://eprints.whiterose.ac.uk/id/eprint/149121/1/Macgregor_et_al_2019_Ecosphere.pdf), [Macgregor author archive](https://zenodo.org/records/1743917), [NERC record](https://nora.nerc.ac.uk/id/eprint/522598/).

Independent review confirmed the outcome distinctions. The148 undamaged successful capsules are across all treatments, not the full-night comparison alone. Code lines49/197 set LitUnlit reference to Lit: a future Unlit coefficient must be sign-reversed for a light/dark direction. Binomial coefficients are log odds ratios, not lnRR. The primary code includes a run intercept but not a plot or reused-individual intercept in these main models; a usage-count covariate does not remove all repeated-individual dependence.

Current boundary: no new pooled pollination estimate, no inferred zero effects for missing extraction, and no wet-versus-dry reproductive moderator estimate. Resolving numeric reproductive contrasts and their dependence is required before that synthesis.
