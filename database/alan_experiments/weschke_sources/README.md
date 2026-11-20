# Weschke2024: prolonged-exposure reef-fish guild richness

[Primary article](https://doi.org/10.1111/gcb.70002), [author data and code](https://doi.org/10.5281/zenodo.14266460). Retrieved and checked2026-11-20. All23 small archive files were acquired and matched published MD5 and byte size; SHA256 fingerprints are also retained. Large example videos were not downloaded or reannotated. Full Methods/Results were read from the cached primary XML, PMC11653166.

## Comparison

Moorea lagoonal reefs, January–April2021. Ten ALAN sites versus ten control sites, pseudo-random allocation with starting-day balancing. Nightly underwater LED exposure was approximately12h; prolonged exposure averaged25 nights. Original camera data are MaxN; following author instructions, each positive taxon entry is converted to presence, then summed within the six author-assigned trophic guilds. This measures detected taxon richness, not fish abundance, feeding activity or predation rate. Some identifications are above species level; no new taxonomy was imposed.

Only16 sites have all three observation periods. The four additional sites (Infront PB, Cooks pharmacy, Club Med Far, Back row purple) have prolonged observations only. These are retained for the final-timepoint comparison, not filled with zero baselines. Repeated-survey and20-site exports use different letter IDs for controls; do not join solely by letter.

## Raw versus derived source discrepancy

Of120 final-site×guild cells,119 match the author-derived guild CSV in the corresponding row order. At Lara Tipanier under prolonged ALAN, the workbook yields2 planktivore labels but the derived CSV gives3. The workbook has positive Dascyllus trimaculatus and Ostorhinchus novemfasciatus values; Amphiprion chrysopterus is0. The separate derived presence table ALAN.Prolonged.B codes Amphiprion as1. This explains the extra species in the derived guild count, but does not establish which source value is correct. Preserve both; do not repair the workbook.

An independent agent reproduced the120-cell audit and localized this taxon discrepancy. It also compared48 repeated surveys for98 exact-matching taxon labels: the same cell differed. Two spelling-mismatched labels were not automatically equated. The reproducible script here audits the20-site guild cells; it does not claim to reproduce that broader98-label check.

## Conditional extraction from workbook

| Guild | Mean ALAN/control | Response ratio | Approximate95% interval |
|---|---|---:|---|
| Piscivores |0.7/0.3|2.333|0.657–8.285|
| Mobile-prey invertivores |2.1/1.2|1.750|0.887–3.453|
| Omnivores |0.5/0.1|5.000|0.464–53.863|
| Planktivores |2.7/1.7|1.588|0.947–2.664|
| Herbivores/detritivores |0.2/0|undefined|not calculated|
| Sessile-prey invertivores |0.2/0|undefined|not calculated|

All guilds are retained. Zero control means receive no pseudocount and no invented finite lnRR. Four finite contrasts use arm-specific sample covariance from ten sites per arm, a first-order delta approximation and Welch–Satterthwaite t intervals. Their full4×4 covariance is retained. Sparse, zero-heavy counts make these approximate intervals uncertain. The author-derived planktivore mean is2.8 rather than2.7; displayed estimates are explicitly workbook-conditional, not a corrected reproduction of the original model.

The point estimates are all positive but individual intervals include1. This does not contradict the original multivariate community test, which has a different estimand and model. Species richness, community composition and predation cannot be substituted for one another. Neither a baseline-adjusted treatment effect nor wet-versus-dry moderation was estimated here.

## Scope and remaining adjudication

Underwater illumination is not a test of water-surface reflection. Reported13–18lux at1m were measured at separate calibration locations, not each experimental site. The study selected anemonefish sites used by parallel ALAN studies; programme overlap must be resolved before cross-paper pooling. Sites were8m or more apart across treatments, so fish movement and shared local environments remain considerations even though light-spillover checks were reported.

Final catalogue additions:0. Four finite conditional contrasts and two undefined contrasts are available for review. Final inclusion requires source-version adjudication and programme matching; these cannot explain terrestrial wet-plant decline directly. Reproduce with `python database/alan_experiments/weschke_sources/extract_prolonged.py`.
