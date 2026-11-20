# Knop2017: model-based light-by-visitor-access interaction

Source: [original article](https://doi.org/10.1038/nature23288), Extended Data Table3 and Methods, 'Estimation of pollination service' and 'Fruit set of C. oleraceum'. Original reproduced article in the government-hosted compilation recorded in source_manifest.json; table visually read from local original_article_complete.pdf page13 (zero-index12), Folio1203. Numeric table cells transcribed, not digitized from a plot. Checked2026-11-20.

## Extractable estimand

The printed interaction treatLamp:baggingOpen is -0.76 (SE0.309) for developed-fruit count and -1.40 (SE0.592) for developed-fruit proportion. These are two dependent endpoints from ONE experiment. They compare the lighting response in open plants against that in bagged plants on the respective model scales. Negative interactions indicate a more negative lighting response with visitor access than under exclusion. They are not the lighting main effect in open plants and are not raw mean ratios.

The fruit-count model is Poisson, and fruit-proportion model binomial; normal Wald intervals are calculated on the printed coefficient scales only. No log-risk ratio or odds ratio is adopted here, as the source Methods does not explicitly name the links. Neither endpoint is pooled with raw-arm lnRR effects. The table includes a lighting main effect for the bagged reference category; calculating the open-plant lighting effect as a coefficient sum requires its covariance with the interaction. That covariance is unavailable. Do not add SEs or assume independence.

## Design and interpretation

In2015, five lit and five dark sites each supported five plant pairs, with bagging randomly assigned within each plant pair before flowering. Total100 plants; two flower heads per plant sampled. Lighting replication is five site pairs, NOT25 independent plants per treatment or200 independent flower heads. White0.5mm mesh excluded flower visitors throughout flowering; it was not a night-only exclusion. The experiment therefore tests moderation by visitor access, not an isolated moth effect, a nocturnal-only contribution, or a formal causal mediation coefficient. A bagging microenvironment effect cannot be separated using this contrast alone.

The original hierarchical GLMM includes plant pairs nested in sites and site pairs, plus an observation-level random effect for overdispersion. Models were selected by AIC-based stepwise deletion; the count model retains herbivore occurrence. Keep model-selection and post-treatment covariate caveats with these published coefficients; no claim of refitting the full model is made. Wald intervals based on rounded published values do not incorporate model-selection uncertainty or small-site-sample correction.

This supplies finite reproductive MODEL coefficients with reported uncertainty, while the existing raw-arm lnRR catalogue remains unchanged. Missing extractable raw-arm effects must no longer be described as a lack of any quantitative reproductive evidence. Between-endpoint covariance and potential programme overlap with Giavi2021 remain unresolved; do not count the two endpoints or companion reports as independent studies. Neither this comparison nor its source establishes wet-habitat moderation or explains the58 candidate species directly.

Regenerate derived intervals with `python database/alan_experiments/knop_sources/compute_reproductive_coefficients.py`.

Independent visual re-extraction confirmed both printed estimates and SEs, reference coding,5site pairs, and all-visitor exclusion. Separate calculation agreed with the coefficient intervals. Standard log/logit links would permit different ratio interpretations, but these assumptions are not needed for the reported coefficient-scale results and are not imposed in the extraction.
