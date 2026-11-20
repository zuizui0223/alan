# Lu 2026: original-data extraction and comparability audit

2026-11-20. Primary result (wet-affinity enrichment) unchanged. This is an independent experimental evidence extraction, not an additional macro result or pooled meta-analysis. No existing meta-analysis effect inputs were used.

## Source and effort verification

[Primary article, doi:10.1002/ece3.74368](https://doi.org/10.1002/ece3.74368); [original supplement package](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles). Author workbook `Date.xlsx`, its `ReadMe`, original `R analysis-ECE.R`, and cached full article XML were inspected. SHA-256 workbook: `ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f`. All source hashes are saved in the audit JSON.

Methods describe five paired lit/control spatial blocks in an alpine meadow, July–August 2024. Day and night each had four occasions of five 20-minute observations per plot (400 planned minutes). Visits required contact with flower reproductive structures. The original R uses the `visitation` sheet for community/guild outcomes, separately for day and night. The `ReadMe` defines `time` as survey round and `repeat` as sampling replicate.

Daytime data independently match the effort description: **20 unique sessions in every block × arm**, four rounds and five replicates; no missing values in required fields. Thus daytime means per original row can be interpreted as visits per 20-minute survey. Five blocks, not 200 daytime rows, are independent replicates. The night rows cannot receive the same clean status: block-arm counts range 18–22, with three duplicated session keys (six rows). They remain flagged pending clarification.

## Recalculated daytime block effects

Calculated from the original community sheet only, with no zero imputation, deduplication, continuity correction, or joining across worksheets. Each block contributes its arm-specific mean of 20 observations. SD and paired arm covariance use the five block means. lnRR is log(mean lit / mean unlit); its delta-method variance includes the paired covariance. These small-sample descriptive estimates are not replacements for the article's fitted-model contrasts.

| Endpoint | Lit mean /20 min | Dark mean /20 min | lnRR | Paired delta variance |
|---|---:|---:|---:|---:|
| Author total column | 16.58 | 17.23 | −0.038455 | 0.007252 |
| Diptera | 7.71 | 5.27 | 0.380488 | 0.006760 |
| Hymenoptera | 5.84 | 7.30 | −0.223144 | 0.005826 |
| Lepidoptera | 0.40 | 0.83 | −0.729961 | 0.290330 |
| Coleoptera | 0.64 | 1.67 | −0.959111 | 0.064809 |
| Other | 2.00 | 2.23 | −0.108854 | 0.099581 |
| Diptera + Hymenoptera + Coleoptera | 14.19 | 14.24 | −0.003517 | 0.005760 |

The three-order aggregate was constructed explicitly by summing those original columns within observation before aggregation. Its mean ratio is nearly one, despite opposing responses among orders. It is not the author total and is not total pollination service.

## Exact unresolved discrepancies

The author total differs from the five guild-column sum in **20 daytime rows** and three night rows. Across daytime observations the author total is 3381 but the five guilds sum to 3389. The workbook Other count is 423, whereas the article reports 482. These inconsistencies have not been silently corrected. Individual guild columns and the explicitly defined three-order sum are reproducibly extractable; treating their sum as a perfectly reconciled community census awaits source clarification. All affected original Excel row numbers are retained in `guild_sum_discrepancies.csv`.

Do not borrow night counts from `species visitation` or `plant–pollinator network`: they are different tables with unreconciled totals and species coverage. Do not treat the companion Lu JPE article as an independent experiment until plots/year overlap is resolved. The original supplement documents observation counts, not visit-specific pollen delivery, stigma pollen, fertilization, or reproductive output.

## Comparison with Giavi 2021

Giavi's [original dataset](https://zenodo.org/records/4540407) records daytime visits to 21 selected plants in three orders, in six lit/dark meadow pairs. Lu's community sheet covers the surveyed flowering community, and its three-order aggregate standardizes order coverage **only**, not plant selection. Giavi uses equal-effort paired sampling events; Lu uses 20-minute quadrat surveys. Within-study ratios remove a common multiplicative time denominator, but do not make plant coverage, spatial sampling or survey windows identical.

Consequently the two studies can be compared as separate daylight visitation responses, preferably by order with sampling design/plant selection documented. They are **not pooled here**. Lu's increase in flies and decreases in other guilds support compositional replacement of visits; neither these counts nor a near-zero aggregate proves compensation of moth-mediated pollen delivery. They also do not establish why wet-affinity plants are overrepresented among macro decline candidates. This experiment is one meadow system, not a wet-versus-dry contrast.

## Reproducible outputs

Under `database/alan_experiments/lu_sources/`:

- `extract_original.py`: reproducible extraction from the cached official source, verifies workbook hash.
- `visitation_original_rows.csv`: unchanged source fields and Excel row identifier.
- `block_totals.csv`: block totals, row counts and planned/implied effort (implied effort is conditional for night).
- `paired_block_descriptive_effects.csv`: block means, SDs, paired sample covariance and lnRR variance; source-row means and block totals are alternative scalings of the same data, **not separate effects**.
- `daytime_effect_covariance.csv`: joint delta covariance of seven daytime endpoints; endpoints overlap and must not be treated as independent.
- `duplicate_session_keys.csv`, `guild_sum_discrepancies.csv`, `extraction_audit.json`: source-integrity barriers and provenance.

Master evidence/effect tables were not modified. Night data and cross-study pooling remain on hold; daytime source extraction is complete with explicit consistency flags.
