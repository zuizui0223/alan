# Deitsch2023: attack proxy and conflicting reported magnitudes

Primary source: [published article](https://doi.org/10.1098/rspb.2023.0153), archived NSF PDF recorded in source_manifest.csv. Rechecked2026-11-20: Methods(Data analysis), Results3a, Figure3 caption and Table1. Page5 visually checked, not solely text extraction. Reported values saved in reported_attack_summaries.csv; no plot-digitized values.

## Reconciliation required

Results3a gives control0.44(SD0.11) and experimental0.53(SD0.09), whereas Figure3 caption gives45.0±10.0% and56.0±9.0%. Differences cannot all be rounding of the same values. Do not average versions or select the larger contrast. The caption does not explicitly identify the dispersion for its arm summaries, so it is not silently labelled SD.

Table1 labels the positive0.09(SE0.03; reported CI0.03–0.16) coefficient 'treatment(control)' but describes increased attacks on experimental plots. The reference coding is unresolved without code/raw data. Preserve the printed coefficient and verbal direction separately. The magnitude0.09 on an untransformed fraction scale corresponds to9percentage points, not a9% relative increase. The text's mean paired ratio1.27(SD0.34) is a third estimand: mean of ratios, not ratio of arm means.

## Independence and endpoint

Each arm's16 samples are four plots observed during four periods. The four plot pairs, not16 plot-period records, provide independent spatial replication. The original attack model pooled spatial groupings after nonsignificant preliminary comparisons and used treatment and sampling period as fixed effects. A nonsignificant spatial test does not prove repeated plot observations independent. Do not obtain a meta-analysis SE by dividing the pooled SD by sqrt(16), and do not convert the selected-model CI into an independent-plot variance.

The endpoint is the fraction of recovered plasticine models with an arthropod mark during exposure. It is not successful killing, adult moth mortality, or a per-hour hazard. Lost/damaged models were excluded by the authors, not counted as attacked. Trapped parasitoids and predators are separate outcomes; their abundances cannot establish successful parasitism or explain the attack contrast by mediation.

## Decision

Retain the reported positive direction as primary-study evidence, with unresolved magnitude and variance. Main effect catalogue remains unchanged. Raw access remains unavailable; the earlier unauthorized EDI API route was not retried. Needed: plot-round recovered and attacked counts, pair/treatment allocation, and a resolution of text/caption/model coding differences. This supersedes any interpretation of the extraction hold as merely a missing-download problem: numerical and dependence adjudication are also required.

Independent visual and Methods review confirmed all conflicting arm summaries, Table1 sign/label, four plots per arm, and the distinction between percentage points and relative change. No independent raw-data reanalysis was possible.
