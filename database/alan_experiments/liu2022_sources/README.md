# Liu 2022: plant growth and insect access

Primary: https://doi.org/10.1111/gcb.16126; original data: https://doi.org/10.5061/dryad.547d7wmb0.
Primary PDF recovered from https://ira.agroscope.ch/de-CH/Page/Einzelpublikation/Download?einzelpublikationId=55545 .

Outdoor Konstanz experiment: 20 cages, 2 lighting treatments × 2 insect-exclusion treatments, five independent cages per cell (Methods 2.2–2.3). Twenty plant species share cages, not 20 independent experiments. Light was approximately 24 lux, cool white 6500K; one dark cage had measurable 0.58 lux and was retained in the main author analysis. Complete versus partial netting manipulated insect access; this does not isolate pollination, herbivory, or predation. Experiment duration was ten weeks; reproductive traits were not measured.

`reported_tests.csv` preserves all 15 fixed-effect likelihood-ratio tests from Table 1 (p3679), square-root transformed aboveground pot biomass. Light × insect exclusion: chi-square=0.009, df=1, P=.927. Commonness × light: 6.779, df=1, P=.009. The author reports biomass responses of +9.2% for common species and −2.0% for rare species (Results p3678). No detectable insect modification is not evidence of equivalence or absence of insect mediation.

The extraction script reads the primary text and records hashes. Independent AI source review confirmed design and key tests; this is not human double extraction. Figure 2 is modeled means ± SE, not raw cage means/SD. No signed coefficient, lnRR, interaction variance, or zero effect was inferred from tests. Dryad metadata identifies the CSV and README but their download endpoints returned HTTP401 on 2026-11-20. Zenodo6007775 contains RMarkdown analysis code only. Raw acquisition remains incomplete.

Unresolved source issues: printed coordinates have invalid sexagesimal minutes/seconds; do not silently normalize them. Reported total 660 pots conflicts with 13 cages × 35 plus seven × 34 (=693). Cage count remains20. Resolve from original CSV before extraction. No effect catalogue rows or current meta estimates changed.

2026-11-21: the public dataset landing page supplied download links `/downloads/file_stream/1343897` (CSV) and `/downloads/file_stream/1343899` (README). Both returned403 through direct retrieval, independently of the previously attempted API401 route. No raw file was recovered. This is an access failure, not evidence that the public files or data do not exist; do not repeat these routes without changed access information.
