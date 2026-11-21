# Frog light × noise experiment: reported model contrasts

2026-11-21. [Primary article](https://doi.org/10.1093/beheco/arac077); [author dataset](https://doi.org/10.5061/dryad.4j0zpc8fj). This extends the existing batch 5 review; it is not a new paper in screening counts.

## Design and distinct response populations

Barro Colorado Island, Panama, July–August 2019: 32 artificial water-filled dishes were installed, 19 became active, and 16 active dishes entered a randomized balanced light × noise crossover in four blocks. Every dish received all four treatments. The same 16 dishes are not 16 independent lit plus 16 independent dark replicates. Light averaged 1.53 lux (range 1.06–3.1); recordings ran 18:00–23:00. The setup selected colonized breeding sites and did not test initial colonization of all 32 dishes.

The paper reports 3,201 recording periods for breeding presence (16 sites, 22 nights), 877 for calling presence, 82 arrival times, 44 calling-start observations (10 sites), 186 calling-effort periods (11 sites), and 122 bouts (8 sites). These are different selected subsets, not interchangeable denominators or independent treatment replicates. Calling effort is the occurrence of a bout of at least 10 calls with no pause exceeding 10 seconds, assessed in selected periods with a single calling male. It is not the number of frogs, total nightly calling, successful mating or recruitment.

## Quantitative results retained

All 28 intercept, light, noise and interaction coefficients from Tables 1–2 are preserved in `reported_coefficients.csv`, along with the original SE, likelihood-ratio test and P-value strings. Positive, null-compatible and negative coefficients are retained. For three explicitly logit-model responses, the reported light coefficient and SE yield the following **single-study model odds ratios**, with approximate unadjusted Wald 95% intervals:

| Response | Odds ratio | Interval |
|---|---:|---:|
| Breeding site presence | 0.835 | 0.381–1.829 |
| Calling presence | 1.062 | 0.128–8.818 |
| Calling effort in selected single-caller periods | 6.172 | 2.362–16.125 |

These use the authors' reduced models after removing the light × noise interaction. They are not unadjusted lit/dark proportions, risk ratios, or three independent studies. Approximate intervals inherit rounding, random-effect specification and model-selection uncertainty; that selection uncertainty is not included in the reconstructed intervals. The author models use date nested within site, plus days active for the first two responses. No covariance between endpoints is available, so the contrasts are not jointly pooled. Different significance outcomes are not a formal test of differences between endpoints.

The arrival-time interaction is reported as −0.65 (SE 0.37; LR χ²=5.10, P=0.02). Methods specify full-model coefficients when interaction is significant, but the table caption broadly calls non-interaction coefficients reduced-model estimates. Arrival-time coefficient provenance therefore needs author-code checking before deriving treatment-specific times. No mean minutes are reconstructed from the rounded intercept and coefficients. LR P-values are preserved, not replaced by Wald P-values.

Only 14 amplexus events and no predation attempts were observed. This does not estimate a null lighting effect on predation, successful reproduction, or population decline. It also does not identify a prey-mediated cause for altered calling. The experimental habitat is water-filled, but there is no dry-habitat treatment.

## Source and acquisition audit

The locally cached full text `PMC9735234.xml` was read again and its SHA-256 is recorded in `audit.json`. Extraction traverses the source tables, preserving their response row spans. Run `python database/alan_experiments/frog_noise2022_sources/extract_tables.py PATH_TO_XML` to regenerate both CSVs.

Dryad's public metadata returned HTTP 200 for version 187596 and lists `activity_data.csv`, `call_measurements.csv`, `call_start.csv`, `site_arrival.csv` and `README.docx`. The activity-data API download returned HTTP 401 on this pass; no access control was bypassed and no raw data were acquired. Metadata is cached externally under `alan_experiment_sources/frog_noise2022`. We have not verified the raw night–dish allocation, equipment-failure repeat, exclusions, or author model execution.

These 28 report rows and three derived contrasts form a separate source table. They add zero entries to the current 108-row group-mean lnRR catalogue and do not change the separately scoped 86-row reported-statistics registry. The evidence is suitable for an explicitly labelled single-study Results subsection, not an added independent predation or pollination meta-analysis.
