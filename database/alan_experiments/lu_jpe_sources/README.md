# Lu JPE2026: recovered publisher supplement

Primary DOI https://doi.org/10.1093/jpe/rtag160 . On2026-11-21 the indexed publisher landing page exposed its own signed DOCX link, which downloaded successfully (198613bytes). The unsigned asset returned403. Full main PDF remains unrecovered. This is supplementary-table review, not full Methods review, final inclusion, or an additional independent experiment confirmed against Lu ECE2026.

`extract_supplement.py` preserves all47 numeric table rows: S2(18), S3(5), S4a(12), S4b(12). Blank repeated-label cells are retained as in the Word table. Strings, direction, SE, df and reported P are preserved, without inventing a sample size or sampling covariance. The table rows are dependent reported statistics, not47 meta effects.

`plot_reported_contrasts.py` forward-fills species/response labels within S4b into a separate derived table, reverses Unlit−Lit to Lit−Unlit, and calculates approximate marginal95% t intervals from rounded SE/df. No multiplicity adjustment, raw-scale back-transformation, or covariance imputation is made. The supplemental figure is a single-report display, not a meta-analysis. Source tables remain unchanged.

2026-11-21 main-source follow-up: Crossref metadata supplied the accepted-manuscript PDF URL https://academic.oup.com/jpe/advance-article-pdf/doi/10.1093/jpe/rtag160/68754751/rtag160.pdf . Direct retrieval returned403; the web reader returned an inaccessible publisher watermark redirect, whose issued link returned400. Crossref metadata is cached externally as `lu_jpe2026/crossref.json`. This establishes the publisher file location, not fulltext access or verification of Methods. Do not repeat this failed route without new access information.

S2 compares pollination access within unlit conditions, not ALAN effects. Its note specifies arcsine square-root seed-set transformation and FDR-adjusted post-hoc P. S4b does not independently state the seed-set scale; application of that transformation to S4 requires the main Methods. Neither table's seed coefficients should be called percentages, percentage-point changes or lnRR.

S4b explicitly uses Unlit minus Lit. For P. kansuensis seed set, Open is0.0445(SE0.0319,df92.8,P.166), Diurnal is−0.0686(SE0.0401,df96.8,P.090), and Nocturnal is0.1332(SE0.0429,df95.5,P.0025). Night-only reproduction differs; all-access reproduction is not clearly reduced. This distinction prevents an abstract-level species-wide decline claim. Lit-minus-Unlit would reverse estimate/t signs only, but we preserve source direction.

S4a reports a light-by-pollination interaction for P. kansuensis pollen(F4.83,P.0099) and seed set(F6.11,P.0031), with no df supplied. E. regelii interaction tests are pollen F.47,P.6270 and seed F.060,P.9414. Differences in within-species significance do not test a between-species interaction. No compensation amount or general self-incompatibility moderator can be estimated from these tables alone.

S3 is visitation per20minutes with plot random effects; its heading says Unlit vs. Lit without an explicit subtraction. Direction remains unharmonized pending main text/code. Its note lists Hymenoptera but the five data rows contain Total,Diptera,Lepidoptera,Coleoptera,Others. Do not infer a missing order to be zero.

Independent AI review read the DOCX tables directly and confirmed the extraction, endpoints and limitations. It is not human independent double extraction. Shared plots/control, independent lighting n, links to the ECE programme, model scale and raw observations remain unverified. No108-catalogue row or pooled estimate changed. Supplement discovery does not increase the full-primary-read count.
