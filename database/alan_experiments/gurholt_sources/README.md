# Gurholt et al. 2026: primary paper and original data recovered, extraction pending

[Paper](https://doi.org/10.1002/ecy.70382), [Cornell original deposit](https://doi.org/10.7298/PARD-TT61). Retrieved2026-11-20. Original294-row data, author README and RMarkdown analysis were obtained through the public institutional repository; deposited byte counts and MD5 checksums match. At initial acquisition the main article and supplement were unread. The dated primary-paper update below supersedes that access hold; original observations and earlier audit findings remain preserved.

The abstract reports faster larval growth and larger adult fresh mass under lighting. Retain this positive-response candidate in screening; effect direction is not an eligibility criterion. No independent effect has yet been estimated.

## Original-data findings

294 rows,96 reused individual labels, unique Cohort+Individual keys. There are30 Cohort+Cage combinations but only8 reused Cage labels; these are not30 independent lighting installations or4 independent studies. First cohort has only Asclepias incarnata (inc), fourth only A.syriaca (syr); cohorts2/3 contain both. Physical cage/light reassignment across years remains to be read.

One row's Treatment disagrees with its Cage prefix; exact ID is retained in `treatment_conflicts.csv`. The author's script derives treatment from the first two characters of Cage rather than from Treatment. Neither is certified as correct without design evidence. Original values remain unchanged.

43 rows lack all six selected fields: initial mass, day7 mass, pupation date, eclosion date, death date and adult wet mass. These are missing observations, not confirmed deaths or zeros. The author's code labels a missing eclosion/pupation date as died; this cannot by itself establish mortality or randomized starting denominators. Missingness must be resolved before survival extraction.

Observed initial/day7/adult fresh/adult dry mass counts by cohort are60/57/24/24,80/46/33/35,80/71/52/0,30/0/13/0. Cohort4 contributes no day7 mass, and cohorts3/4 no dry mass. Thus the same four cohorts cannot support every endpoint. Both adult-wet-mass columns coincide at a factor1000 on121 complete pairs; do not concatenate or average their numeric values as if measured in the same units. Documentation/code must establish units.

The analysis code creates CohortCage random intercepts and growth.rate=(Day7Mass-InitialMass)/6. This code inspection is not model reproduction or validation of the six-day denominator for every row. It does not settle dependence across repeated physical cages. Code package installations were not executed.

## Extraction decision

Hold effect extraction until full lighting design, treatment mismatch and missing-data meanings are resolved. Preserve host-specific support and stage-specific endpoints; fresh mass, dry mass, wing area, developmental duration and survival are not interchangeable outcomes. Host-plant rearing is not a measurement of pollination or predation. No main effect-catalogue rows, new pooled estimate or primary-read count added.

Run `python database/alan_experiments/gurholt_sources/audit_raw.py` to regenerate support counts and conflict records. Source hashes and original-column names are preserved. This is acquisition and data-structure verification, not completion of eligibility review.

## Independent data/code check

A second reader independently confirmed the support counts and identified the conflict as cohort3 IndividualML4L (TreatmentAL, CageML4). Initial+day7 mass coexist on174 rows, all cohorts1–3. Code converts Adult wet mass to mg by multiplying1000, supporting g-versus-mg interpretation, but the two mass columns have different coverage (122 versus154 observations;121 overlap), so no missing values were filled across them. README lists2024 logger prefixesWL while raw individual Cage labels areAL/ML; their physical crosswalk remains unverified. Eclosion date is present on130 rows; of164 without it,106 also lack death date. Death annotations include missing individuals, and four rows contain both eclosion and death dates; simply coding no eclosion as larval death would misclassify unknown/other-stage records. These checks reinforce the extraction hold rather than licensing a survival estimate.

## Primary-paper recovery and design adjudication, 2026-11-20

The13-page published PDF was recovered from the [corresponding research laboratory publication list](https://agrawal.eeb.cornell.edu/publications-1/). Methods/Results pp3–6 read. Supplement remains unrecovered.

**Confirmed field design:** cages along a transect9.14m apart, every other cage assigned a white solar LED (~7lux,178cm aboveground and30cm from cage); controls natural night (~0lux). Assignment is alternating, not stated as independently randomized. Location changed from42.503274,-76.467617 in2023 to42.443962,-76.501884 in2024. Cohort1 used3lit/3control cages; cohorts2–4 used4/4. These cohorts are repeated within year; do not count30 cohort-cage observations as30 independent installed lamps. Cage/year is the candidate illumination unit and repeated2023 cohorts need shared-unit covariance. A physical ID diagram and exact reuse/relocation within2023 remain to be checked, rather than assuming30 independent units.

**Starting support:** Methods states10 larvae/cage in2023, hence60/80/80 larvae in cohorts1/2/3, consistent with initial-mass rows but not the72/96/96 template rows. Cohort4 used64 plants but30 larvae (15per treatment), so plants must not be substituted as larval denominators. This helps explain surplus template rows without assigning individual blank rows a biological fate. Thirty2024 larvae matches the raw table; paper says14adults at completion, versus13adult-wet-mass values in that cohort.

**Endpoint and version differences:** paper reports127fresh masses and62dry masses, while deposited columns have122 and59 observed values. Alternative adult wet mass.1 has154 values and different coverage; do not silently replace the first column to force agreement. Paper describes a7day growth denominator, while code uses6 and174complete initial/day7 pairs. The common fixed divisor would cancel in an unadjusted treatment mean ratio, but does not resolve sampling or individual intervals, nor authorize an invented correction.

**Reported findings, not independently reproduced effects:** Results text gives growth +16% with p=.053, larval duration -3% with p=.128, pupal duration +5% with p=.040, and fresh mass >7% higher with p=.010. Dry mass p=.315 and wing loading p=.200 do not establish increases. Figure1 caption instead lists growth p=.038 and pupal duration p=.029; preserve rather than choose the more significant value. Path-analysis growth p=.028 is a distinct model, not resolution of the text/caption discrepancy. Host-by-light wing-area interaction is reported p=.035. Neither survival/reproduction benefit nor pollination is established by these growth and morphology results.

The abstract's wording about consistently faster growth is therefore insufficient for a significance claim across endpoints/models. Register as a primary-read field experiment with raw data; keep effect extraction pending explicit endpoint/cohort selection, allocation discrepancy and sample-count reconciliation. Do not exclude because some directions are positive, and do not interpret higher fresh mass without higher dry mass as enhanced biomass production or population growth.

Independent primary-paper review confirmed the design and counts above. It also found Fig.2 fresh-mass p=.002 versus Results text .010, and host-by-light wing-area p=.031 versus text .035. Initial-stage/move dates in the raw data include second-instar/next-day records despite the Methods description of first-instar/same-day transfer. These are additional reconciliation items, not silently corrected data.

## Growth interval adjudication, 2026-11-20

`audit_growth_intervals.py` verifies the unchanged deposited file checksum and exports all174 complete initial/day7 mass rows and22 cohort-cage means. All174 recorded day7 dates are exactly7days after the move date. Intervals from hatch are7days (36),8days (136),10days (1), and -357days (1). The negative interval is cohort1 ML2H: hatch2023-06-22 but move2022-06-23 and measurement2022-06-30. Original dates remain unchanged; the correct year cannot be inferred from arithmetic alone. Initial weighing has no separately recorded date, so this confirms the move-to-measurement interval, not the actual initial-weighing interval.

Normalized move stages are first instar119, second instar54, unknown1. All mass increments are positive. The previously identified cohort3 ML4L treatment conflict has neither mass and is absent from the growth subset; all174 complete growth rows agree between Treatment and Cage prefix. Therefore that specific conflict is NOT a reason to block this endpoint, although it remains relevant to other endpoints.

On identical cohort-cage aggregation, dividing every mass increment by6 or7 yields identical log mean ratios and delta-method variances (maximum numerical difference2.1e-16). This is an algebra verification, not an adopted effect or validation of sampling independence. Do not divide by hatch-to-measurement intervals: the initial measurement was not established as hatch mass, and one hatch interval is invalid. The remaining growth extraction questions are endpoint/cohort selection, stage composition, initial weighing time, and dependence of reused physical cages. The22 cohort-cage means are not22 independently installed lamps. No effect was added to the main catalogue or pooled model.

An independent reader separately reproduced the174 complete rows, both interval distributions, stage counts, absence of treatment conflicts, and positivity of increments, and agreed with the common-divisor invariance argument. Reproduction: `python database/alan_experiments/gurholt_sources/audit_growth_intervals.py`. Row-level provenance is in `growth_interval_rows.csv`; summary and algebra checks are in `growth_interval_audit.json`.

## Descriptive cohort mass-gain contrasts

`extract_cohort_mass_gain.py` now extracts the observed mass increment (follow-up minus initial) among complete records. It uses equal weights for cage means, no day divisor, and all three supported cohorts. It does not correct original dates, substitute adult columns, infer survival or resolve host/stage composition. Cohort4 has no complete pair and is not assigned a zero effect.

| Cohort | Complete individuals | Cages per arm | Lit/dark ratio | Approximate95% interval |
|---|---:|---:|---:|---:|
|1|57|3|1.073|0.975–1.181|
|2|46|4|1.321|1.015–1.718|
|3|71|4|1.036|0.791–1.357|

Intervals use between-cage delta variance and Welch-type degrees of freedom within each cohort; no multiplicity correction. They are not the author's adjusted model intervals. Cohort2's interval excluding1 is not an independent replication of a general lighting benefit. Cohort1 is mostly second instar at transfer; cohorts2/3 mix host species. Selection into two observed weighings and host/stage imbalance remain possible explanations of the descriptive contrast.

All three cohorts were conducted in2023 and reuse cage labels. Their cross-cohort covariance is unknown, not zero; no pooled estimate is computed. These are three dependent comparisons from one programme, not three new studies. The main catalogue remains unchanged. Independent recalculation reproduced all174 complete rows, per-arm cage/sample support and all three mean ratios from the original file. Outputs retain cage means and source hash for audit.
