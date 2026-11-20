# Gurholt et al. 2026: original data recovered, design pending

[Paper](https://doi.org/10.1002/ecy.70382), [Cornell original deposit](https://doi.org/10.7298/PARD-TT61). Retrieved2026-11-20. Original294-row data, author README and RMarkdown analysis were obtained through the public institutional repository; deposited byte counts and MD5 checksums match. Main article Methods/Results and supplement remain unread. Do not advance this report to fully reviewed on the strength of the abstract or original files alone.

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
