# Independent exploratory leaf-damage synthesis review

2026-11-20. Reviewed the pre-computation specification `leaf_damage_meta_plan.md`, existing source-design notes, original Hong Kong herbivory CSVs, and Cieraad/Ceppi extracted contrast/covariance matrices. This independently verifies synthesis arithmetic; it is not a fresh complete re-extraction of the Dutch and Swiss raw workbooks or rerun of published models.

Sources: [Guo2025 primary](https://doi.org/10.1007/s00114-025-02031-w), [original Figshare](https://doi.org/10.6084/m9.figshare.27273918.v1); [Cieraad2023 primary](https://doi.org/10.1111/1365-2664.14336), [original Dataverse](https://doi.org/10.34894/VD98DG); [Ceppi2025 primary](https://doi.org/10.1111/1365-2664.70040), [original EnviDat](https://doi.org/10.16904/envidat.628).

## Independent numerical calculation

For Cieraad, included every six-row contrast and the matching 6×6 covariance; for Ceppi, selected exactly the three `all_leaf_mean_percent` rows and their covariance submatrix. Calculated equal-weight mean log response ratios and `a'Va`. For Hong Kong, read `Herb_Psyc.csv` and `Herb_Aporosa.csv` directly, restricted to Sampling=3, averaged Percentage by Site/Setup, then formed a four-variable covariance matrix across the four sites. Applied the log-ratio derivative matrix to that covariance divided by four, retaining covariance between species. Thus the HK calculation did not reuse the prior seven-pair effect file.

| Programme | Mean log ratio | Sampling variance | Exp(mean log ratio) |
|---|---:|---:|---:|
| Cieraad | 0.646205197509722 | 0.123289277256756 | 1.90828550487349 |
| Ceppi | −0.386332519183647 | 0.108049293005311 | 0.679544526516665 |
| Guo, four-site weighting | −0.329085158552935 | 0.006634547693072 | 0.7195817 |

HK arm means are Psychotria lit 4.415, unlit 7.6402125; Aporosa lit 2.874625, unlit 3.2080875. Species log ratios are −0.54841758 and −0.10975273. Their covariance is approximately [[0.16967158, −0.18625104], [−0.18625104, 0.22936869]]. The small variance of their average results from strong negative covariance, not large replication. It must not be interpreted as unusually secure evidence: only four site vectors estimate that covariance.

Independently optimized the scalar REML objective `0.5*(sum(log(v+tau2))+log(sum(w))+sum(w*(y-mu)^2))`, with `w=1/(v+tau2)`, comparing the zero boundary. Obtained:

- tau² = 0.209797747314
- pooled log ratio = −0.073909748354
- Q/(k−1) = 1.082395962250
- modified-KH variance = 0.100512512082
- ratio = **0.928755511773**, t2 95% interval **0.237399662249–3.633479477082**
- normalized programme weights = 0.27878939, 0.29215666, 0.42905394

These numbers are standalone checks for the forthcoming implementation, not validation of files that did not yet exist at calculation time. No normal approximation replaced the t2 interval; no off-diagonal covariance was discarded.

## Estimand and comparability assessment

**Pass for the stated exploratory fixed-menu synthesis, with substantial limitations.** All included endpoints describe retained-leaf area damage, including undamaged leaves, rather than damaged-leaf incidence or severity conditioned on damage. Percent versus fractional units cancel in ratios. Nevertheless, the measurements are not necessarily total consumed leaf area divided by total available area: leaf-level averaging and plant/installation weighting differ. The result is a mean of programme-specific log responses, not a universal tissue-consumption rate.

Each programme averages a different menu of species, spectra, doses, durations and settings. Keeping all six Dutch, three Swiss and two HK contrasts prevents selecting only favourable responses, but does not standardize exposure or species. The equal average of log ratios is a geometric-mean response, distinct from a ratio of pooled arithmetic counts/areas. One programme remains one meta-analytic unit regardless of menu size.

Only two forest/site blocks support each Dutch/Swiss covariance, with rank-one contrast covariance; the four HK sites likewise provide weak variance estimation. REML and modified KH account for neither full uncertainty in these estimated covariance matrices nor bias from retention of measurable leaves/plants. The broad interval covers substantial decreases and increases. Absence of other Dutch/Swiss reports prevents within-this-synthesis duplicate programmes, but does not finish the overall provenance adjudication or systematic search.

Final-round HK selection is fixed rather than chosen for direction; its earlier-season signal differs. Keep Aporosa at genus level pending conflicting source-name resolution. Do not convert no records into zero damage, generalize to unmeasured species, or claim this establishes moisture modification, spider predation, chemical mediation, pollination decline or the cause of the macro-level wet-affinity association. Final systematic eligibility remains false; the result can enter a clearly labelled exploratory Results subsection without becoming the main wet-versus-dry test.

## Independent wet-cohort and trait-scope check

Re-read `EIV_Europe.csv` with Latin-1 encoding and species index, and `results/all_species_forest/wet_reproduction/species.csv`. The latter contains 363 labels: 58 candidate and 305 comparison rows, representing 361 distinct accepted names. Exact-binomial results:

| Experimental species | EIV MOIST | Exact source-name matches | Exact accepted-name matches |
|---|---|---:|---:|
| Fagus sylvatica | 5.5 | 0 | 0 |
| Quercus robur | 5.7 | 0 | 0 |
| Vaccinium myrtillus | x: nonnumeric | 0 | 0 |
| Psychotria asiatica | No exact EIV entry | 0 | 0 |

Aporosa remains unresolved at species level because of the source-name conflict. Its trait/cohort match is **unassessed**, not a confirmed absence. No synonym resolution was performed. The Vaccinium `x` and missing Psychotria entry must not be encoded as dry habitat affinity or numeric zero. Fagus and Quercus numeric moisture affinities do not identify experimental soil moisture or proximity to water.

The leaf-damage synthesis therefore contains no exact-name overlap among its four resolved species with this specific wet-candidate/comparison cohort, and it does not compare high versus low moisture affinity or wet versus dry experimental settings. Its imprecise pooled response is neither a direct test nor a causal refutation of the proposed mechanism behind the wet-plant macroassociation. It provides a separate assessment of recorded leaf damage under the selected illumination menus. The unresolved Aporosa identification and incomplete trait coverage further limit any transfer to that macrocohort.
