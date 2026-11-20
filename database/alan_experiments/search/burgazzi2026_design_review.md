# Burgazzi2026 design review

Reviewed 2026-11-20. Methods and Results read in the [publisher full text](https://nsojournals.onlinelibrary.wiley.com/doi/10.1002/ecog.07758?af=R); primary DOI [10.1002/ecog.07758](https://doi.org/10.1002/ecog.07758).

## Primary-study summary (under 200 words)

In August 2021, 18 branching flumes were arranged in six header-tank blocks. Each block contained control, separated-stressor and co-occurring-stressor flumes. Each flume had two upstream arms joining downstream. Flow cessation and 10-lux, 4000-K lighting were assigned upstream. Drift was sampled over 24 hours weekly; benthos was sampled twice per section at day 22. The Methods describe header/flume nesting; the recovered code differs (see below). Downstream drift was computed by adding upstream catches, not independently sampled downstream.

Light alone increased drift richness but did not significantly change total drift abundance; benthic richness decreased. Flow dominated several responses, with metric-specific interactions. Downstream drift export differed among stressor arrangements, whereas overall downstream benthic richness/abundance did not significantly differ. Predation, emergence and consumption were not measured; predator avoidance and food limitation are proposed explanations.

## Data access and extractability

[Dryad dataset and README](https://doi.org/10.5061/dryad.rfj6q57k5) list `DISCONNECTED_alldata.xlsx` and describe drift, benthos and environmental sheets. Identifiers include header tank (`HT1–HT6`), flume (`F1–F18`), upstream side, treatment and drift sampling date/time. Two benthos samples per section are described. The linked [author code repository](https://github.com/GSBurg/Disconnected) contains `Disc_macroR1.R`. Workbook values were not recovered or audited in this review. The parent source-access audit reports public download 403 and API 401. The recovered author code was subsequently inspected as described below; listed workbook files do not establish successful numeric extraction.

## Own design audit and extraction implications

- Header tank is a shared block; tributaries share a flume. Retain both structures. The allocation implies six light-only, six flow-only and six combined upstream arms, plus 18 untreated arms: two in each control flume and one in each combined flume. Confirm this against raw IDs before extraction. Controls cannot be treated as 18 wholly independent replicates against six treated arms.
- A transparent primary candidate is the light-only upstream response versus an explicitly selected, block-matched untreated comparator, retaining six header-block contrasts. A flow-by-light modification requires the full four-cell covariance, not independent lnRRs. Controls in combined-treatment flumes and in wholly untreated flumes need an explicit averaging rule because their connected contexts differ.
- Collapse or model benthos subsamples within section. For drift, T4 is a possible terminal observation under a last-time rule; retaining all four weeks requires repeated-measure covariance. Do not select the most responsive week or taxon.
- The sum of upstream nets is an export estimand. It cannot serve as independent confirmation of the upstream response, a downstream retention rate, or an abundance census. Benthos standing stock and drift flux must remain separate endpoints.
- Light-only local contrasts should not be relabelled as the downstream comparison between combined and separated stressors: the latter estimates spatial arrangement of two manipulations.
- Suitable in principle for aquatic movement/standing-stock stages after unit-level data verification. It is not a direct predation or plant-reproduction experiment, nor an experimental wet-versus-dry comparison. No effect-size number extracted here.

## Possible overlap with Juvigny-Khenafou2024

The [2024 primary article](https://doi.org/10.1016/j.scitotenv.2023.168836), examined through indexed publisher content and [PubMed](https://pubmed.ncbi.nlm.nih.gov/38016568/), describes 18 Sulzbach-fed flumes at the same research station and flow/ALAN manipulations, measuring litter decomposition and associated invertebrates. The indexed report states limited ALAN effects. Full publisher access returned 403 during this review; exact experimental dates, allocation IDs and raw-file overlap were not verified.

This is a concrete shared-programme flag, not proof of duplicate measurements. Until date/ID matching resolves it, do not count the two publications as independent experiments. If the same campaign is confirmed, retain distinct endpoints within one programme and account for shared controls; apply the database first-report rule to genuinely repeated endpoints. Do not discard the 2024 null evidence because later endpoints respond.

## Recovered code discrepancy (independent inspection)

Inspected externally cached `Disc_macroR1.R`. Lines 40 and 591 build `nesting` by concatenating HT and flume. The univariate models use only `(1|nesting)` (e.g. lines 58, 94, 621, 657), not separate header and flume-within-header intercepts. Concatenation identifies a flume-within-header group but does **not** estimate covariance between different flumes supplied by the same header. It is therefore not equivalent to `(1|HT/flume)`. This differs from a literal two-level interpretation of the Methods description. Original inferential uncertainty should not be assumed to include the shared-header component without clarification.

The upstream abundance models at lines 94 and 657 use square-root-transformed responses. Their coefficients (including Table 2 abundance entries) are not lnRRs; exponentiation would be invalid. Downstream benthos uses `log1p`, also not ordinary log-mean response ratios. Raw observations are still required for the proposed self-extracted estimand. No original model was rerun or repaired here.
