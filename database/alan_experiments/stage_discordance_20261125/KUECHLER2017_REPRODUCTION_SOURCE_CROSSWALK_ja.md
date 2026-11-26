# Newly recovered 2016 Swiss light-by-diel-access experiment: Küchler 2017 (2026-11-26)

## Evidence scope

New grey-literature report found by targeted search, not a new confirmed independent lighting programme. Olivia Küchler, 2017, University of Bern master's thesis, *Effects of artificial light at night on plant reproductive success*. Primary text: https://www.parcs.ch/frg/pdf_public/2022/51204_20220906_150635_2017_MasterThesis_OKuechler_final.pdf. Original text and appendix table text inspected; requested PDF screenshots failed cache miss. No visual double extraction, source raw observations or peer-reviewed counterpart obtained. No digitized figure means. Citation and discovery retained even where results are null or opposite to hypothesis. Conference abstract by Lofchick et al. (2023; https://entomology.org.il/wp-content/uploads/2023/10/Israel-Entomological-Society-2023-abstracts-1.pdf) is another preliminary reproduction lead, not numerical data or a verified independent programme. Tel Aviv University project page https://lifesci.tau.ac.il/botanical/research/Silenepalaestina01 confirms topic; no full experimental dataset acquired.

## Design and independence

Three lit and three dark Swiss ruderal sites (Table 1), two blocks/site, four access conditions (open, day-only, night-only, permanent exclusion), four plant species. 169 plants, 582 inflorescences, 9651 seeds (Appendix 2). Plants randomized to access treatment; physical lamp treatment randomization not established here. Study date text gives 16 July–15 August 2016 but calls it 40 days; calendar interval differs. Do not silently resolve this. Potted C. jacea, D. carota and S. vulgaris versus in-situ C. oleraceum; D. carota and S. vulgaris insecticide exposure, S. vulgaris cutting. Mesh and handling may affect plants. Plant/cage/site nested models, GLMM count model with observation-level random effect; C. oleraceum mean per head treated by LMM. Link and exact model coding unverified.

Geodesic nearest-site match (Haversine, radius 6371008.8 m) against frozen Knop2015 site coordinates (14) and Giavi2016 site coordinates (12). Same light arm required for match; all six thesis sites match both sources within 50 m. Knop distances 5.25–21.75 m, Giavi 4.68–35.14 m. Knop2015 four of these sites included the reproductive experiment, all six included the visitation experiment. Giavi2016 complete paired IDs 1,2,4 coincide, with overlapping sampling period. These are approximate coordinate matches, not proof of identical plot borders or plant identity. Original site locations are rounded and lighting history unclear. Conservatively cluster these reports in the shared Bernese meadow programme; never count six sites, four species or three publications as independent programmes. New report count +1, independent programme count +0.

## All reported seed-number tests

|Species|light main p|access p|light × access p|
|---|---:|---:|---:|
|Centaurea jacea|.0661|<.001|.0860|
|Cirsium oleraceum|.5270|<.001|.7578|
|Daucus carota|.8873|.5467|.6168|
|Silene vulgaris|.6435|.4266|1|

The S. vulgaris interaction has suspicious delta-AIC -59.59 and LRT=0; convergence/model selection need confirmation, not equivalence. All four mean seed mass main and interaction p-values retained in JSON. For C. jacea Appendix 3b, conditional on inferred usual treatment coding, open light coefficient +3.783 ±1.193, day light coefficient +.315 (joint SE range .403–2.789), open-minus-day light contrast +3.468 ±1.596. The nominal normal reference interval for the latter is [.340,6.596]; a t4 reference with six physical sites is [-.963,7.899]. The latter is a sensitivity, not recovered model df. The omnibus interaction is .0860, and the text reports open-light z=3.012 whereas rounded coefficient/SE gives 3.171; source inconsistency retained. No exponentiation or rate interpretation before link/code confirmation. Dark night-only arm has zero seed production and unstable coefficients (SE ~461); permanently caged comparison similarly unstable (SE ~515). Do not impute a pseudocount or report a very large finite ratio. C. jacea positive open outcome and night-only seed presence are informative direction/context; no direct proof of the proposed diurnal-pollinator extension mechanism.

C. oleraceum Appendix 4a: under the inferred reduced access coding, intercept 9.2327 and contrasts permanent -5.3159 ±1.0338, day -1.7417 ±1.0332, night -2.0027 ±1.0332. Model-predicted access-specific values O 9.233, B 3.917, D 7.491, N 7.230. Relative to model intercept, losses .576, .189, .217 respectively. These are model-derived access contrasts across unspecified light reference/averaging, not dark-only measured dependence or light effects. Day-minus-night coefficient .261; unknown covariance gives normal outer [-3.789,4.311]. No certain difference between temporal routes. C. oleraceum is an exact match to frozen wet decline candidate VJ8V, EIV MOIST 7.2. Knop2017 reports decreased fruit set in 2015, whereas thesis does not resolve a seed number effect in 2016; different endpoint and overlapping sites, not independent replication or effect reversal. Thesis discussion itself describes lower-direction tendencies and limited replication.

S. vulgaris flower parasitism: overall light coefficient -.1516 ±.0805 (reported p=.0884); t10 reference [-.331,.028]. Open comparison reported t=-2.143, p=.0386; interaction p=.4734. Different host and Hadena species (H. rivularis), compared with S. latifolia/H. bicruris in Boom/Giavi. Opposite directions of seed-predator effects matter; no shared outcome-scale quantitative pool. No causal net seed benefit demonstrated.

## Co-located daytime visit evidence

Using existing author-derived Giavi2021 original RData paired arm means at matching three sites, equal weight across pairs, compute lit/dark ratio of arithmetic means, paired delta scores and nominal t2 interval. This is community visitation to selected 21 plants, not the same plants or seed observations. Full six pairs included as reference, no independent additional k.

|Response|Full 6 pairs RR|Overlap 3 pairs RR (95% nominal t2)|
|---|---:|---:|
|3-order total|1.052|1.059 [.651,1.723]|
|Diptera|1.004|.918 [.072,11.690]|
|Hymenoptera|1.072|1.285 [.219,7.560]|
|Coleoptera|1.185|.880 [.136,5.683]|
|Diptera/Hymenoptera response ratio|.937|.714 [.012,41.874]|

All intervals include 1, leave-one-pair values saved. No statistically resolved day-guild compensation, no contemporaneous individual mediator data.

## Cross-programme interpretation and admission

Lu JPE S4b open-minus-day light contrast (lit−dark convention): P. kansuensis -.1131 (max-SE normal outer [-.2542,.0280]), E. regelii +.0156 ([-.1100,.1412]). Source scale not harmonized, no pooled effect with C. jacea +3.468. New thesis offers a particularly relevant positive open response rather than a universal loss. The same C. jacea in Grognuz 2026 had seed RR .822 overall, .302 without slugs, 1.289 with slugs (provisional table transcription; reused illumination locations and missing outcomes). Different exposure and seed denominators; do not infer a verified same-species interaction across studies or pool them. Additional field phenotype and seed timing studies without measured reproductive output are separate.

Conclusion: **do not admit a new independent reproductive meta effect**. Add one grey-literature report to intake; update programme crosswalk. Retain access, seed-number, parasitism and mass as distinct endpoints. The new evidence gives a same-year shared-site bridge and an actual wet candidate experiment, but the conditional model covariance and comparable raw outcomes remain missing. Macro decline candidates remain occurrence-rate patterns, not proved demographic loss or causal mediation.

Output kuechler2017_reproduction_source_crosswalk_20261126.json (blob 828e316ecd40eda6f9fb5fa9e568b7127f1985ef), script scripts/analysis/alan_kuechler2017_reproduction_source_crosswalk_20261126.mjs. Six supporting Git blobs pinned, 15 numerical identities checked from saved output (max error 5.551115123125783e-17), all 12 proximity matches verified. Independent CI, primary table images, original plot/plant data, full-source human adjudication not obtained. Main and existing formal catalogue/pooled results unchanged.
