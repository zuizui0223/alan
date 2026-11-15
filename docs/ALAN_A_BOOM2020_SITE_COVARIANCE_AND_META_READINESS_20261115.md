# ALAN reproductive meta-analysis: original Boom 2020 covariance and source-scale sufficiency audit (2026-11-15)

## Original experiment and source

Boom et al. (2020), DOI https://doi.org/10.1038/s41598-020-75471-1, **four original lit/dark spectrum field sites**, four source treatments Dark/White/Green/Red. The original Europe PMC JATS XML Table 1 (PMCID PMC7591485) was re-fetched and its DOI, four source stages, original per-arm printed seed means and **per-arm standard errors (SE)** validated against exact source transcription. Source XML SHA256 `182dafe1e99fef5453ecf0d1d1401361fe949fd8d8617ae040588fd19795154f`. [Original original-source replay GitHub Actions #43724951283 SUCCESS, 5/5 tests](https://github.com/zuizui0223/alan/actions/runs/43724951283). Source script [here](../scripts/analysis/alan_A_Boom2020_original_SE_covariance_20261115.py) and machine-readable numerical report [here](../results/alan_A_Boom2020_source_signed_covariance_meta_readiness_20261115.json).

**Important endpoint**: “estimated seeds per flower” comes from the author's *product of treatment-specific mean fruit set × (1−Hadena larval fruit infestation) × fertilized ovules among selected uninfested fruits × mean 404 ovules*, not total actual mature seeds from every originally marked plant. It is not true total plant seed yield, and conditional ovule estimates select on post-treatment fruit infestation.

| Original spectral treatment | Published estimated seeds/flower ± arm SE | Printed treatment / dark ratio | Possible SE of arm difference with unknown rho [-1,1] | Maximum hypothetical |t| if 4 true independent paired sites |
|---|---:|---:|---:|---:|
| Dark | 117 ± 40 | 1 | — | — |
| White | 86 ± 28 | 0.735 | 12–68 | 2.58 |
| Green | 61 ± 20 | 0.521 | 20–60 | 2.80 |
| Red | 90 ± 28 | 0.769 | 12–68 | 2.25 |

## New quantitative inference constraint

For any two original estimates X (lit) and C (dark), their difference has sampling uncertainty

\[
\mathrm{SE}(X-C)=\sqrt{\mathrm{SE}_X^2+\mathrm{SE}_C^2-2\rho\,\mathrm{SE}_X\mathrm{SE}_C}.
\]

Here \(\rho\) is the missing paired **site-level sampling error correlation**; the published Table 1 reports individual arm SEs but not their shared-site covariance. Even assuming the three spectra all share the same dark reference, the observed mean differences **cannot be given one valid unique contrast SE or meta-analysis weight**. Merely treating all 600–700 flowers as independent assignments would be pseudoreplication.

**Illustrative, NOT the author's fitted model**: if the four field sites were four independent, normally distributed paired units and the reported per-arm SEs were valid site-level SEs for those exact paired units, the absolute paired t statistic could be at most 2.58 (white), 2.80 (green), and 2.25 (red) for any admissible \(\rho\). For df=3, a two-sided 0.05 t threshold is 3.182, so *even the most favorable mathematical correlation* would not pass that conditional 4-site t test. This demonstration is **not a new p-value**, cannot replace the published model or recover the missing covariance, and should not be interpreted as a fully reconstructed cluster analysis. The authors' **published omnibus estimated-seed-output ANOVA p=0.63** remains the actual reported finding.

Stage decomposition from the printed rounded fruit/infestation/ovule percentages yields 0.727 (white), 0.515 (green), 0.760 (red) relative to dark; published estimated seed means yield 0.735, 0.521, 0.769. This discrepancy is **only source table rounding** and should never be promoted to two independent reproductive estimators.

## Meta-analysis implication

This publication is **one** true four-field experiment with **three correlated lighting spectrum contrasts**. It provides (i) stage-specific directional information, (ii) published arm SEs with unknown covariance, and (iii) direct evidence that adult moth pollination and larval fruit predation can respond differently. It does **not** add three independent seed-outcome effects, confirm a significant green-light net seed loss, or estimate a night-dependent ALAN × V meta regression.

Under the strict proposed primary meta-analytic endpoint **open access mature seeds per originally assigned flower and its true light-site-clustered log RR variance**, this paper remains **NOT READY**. It cannot be fixed by inverse-variance pooling the published three arm SEs; the original site-level estimated seed-output pairs, original marked-flower seed counts and shared-site covariance need to be recovered.

The literature audit after targeted updated searches has not established another independent paper with identical total initial-flower seeded-output and correctly clustered variance; the authenticated existing nine direct ALAN reproductive publications remain nine, the observational three remain separate, and the 2023 Silene palaestina conference abstract remains grey literature. No fabricated null (unknown point effect) or invented precision is inserted. A systematic PRISMA search remains to be completed before claiming global exhaustiveness.

This is a post hoc original-source methodological audit and provenance-backed inferential gate, not a registered meta-analysis, not a new causal estimate and not a global empirical conclusion about islands.
