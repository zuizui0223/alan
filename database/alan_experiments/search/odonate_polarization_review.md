# Odonate polarization experiment: daylight scope and version reconciliation

Reviewed 2026-11-20. Both versions' Methods and Results were read directly. **Exclude from ALAN manipulation effect catalogue: no artificial nighttime illumination was manipulated.** Retain as separate mechanistic evidence on surface cues, not as an ALAN effect or a wet-versus-dry experiment.

## Primary sources

- Preprint: https://doi.org/10.21203/rs.3.rs-4458491/v1 ; Research Square PDF downloaded from https://assets-eu.researchsquare.com/files/rs-4458491/v1/39ad595f-71a0-4a74-a78e-14b482419a0c.pdf?c=1720483534 . Posted 2024-06-03, 23 pages.
- Final: https://doi.org/10.1093/beheco/arag108 ; PMC13595255. Final Methods/Results obtained as NCBI BioC XML through its public full-text service. Public PMC HTML returned a browser check; this was not counted as full text. The downloaded XML actually contains Methods, Results and tables/captions and was parsed for reading.
- The same title, ordered five authors, four Mexican localities, 20 choice sites, four treatment surfaces, 136 observed individuals (130 male, 6 female), 17 species and matching temperature means support manual identification as versions of one experiment. Both DOIs already occur in the search inventory. Linking versions does not imply unchanged models or results.
- Final data citation: https://doi.org/10.6084/m9.figshare.33065717 . An earlier search lead points to 27633639; neither dataset has been downloaded in this review, and no identity between those deposits is assumed.

Cached primary PDF/XML/text remain outside Git under `../gbif_annual/alan_experiment_sources/odonate_polarization/`; hashes in the source manifest. No raw observations or new effect estimates were extracted.

## Design relevant to eligibility

At four localities near Chamela-Cuixmala, Mexico, five sites per locality received four adjacent tray treatments: green or stainless-steel surface, with or without water. Trays were placed approximately 50 cm from water, 25 cm apart; sites were separated by 100 m. Positions were randomized. Plots were deliberately exposed to **natural sunlight**, and two 50-minute observations/day were performed between **11:00 and 14:00**. The manipulation changes surface material/polarization and water, not illumination versus dark controls. Surface type also changes appearance and temperature; it is not an isolated manipulation of polarization alone.

Urbanization/polarization-exposure categories were assigned using vegetation and satellite-image car counts (0,10,137,257). These are not measured ambient nighttime irradiance, VIIRS radiance or field measurements of environmental polarization. One locality represents each exposure category, so category is confounded with locality. The four focal species did not all occur in every locality. Daytime preference cannot establish inherited adaptation to chronic ALAN.

The observed responses were perching, patrolling and aggressive interactions, predominantly in males. Individuals were unmarked. Repeated behaviors and the four adjacent trays must not be treated as independent lighting replicates. A female was observed laying eggs in both water treatments; reproductive success was not measured. Dead beetles on the low-polarizing metal surface were descriptive observations, not an ALAN mortality or predator-effect contrast.

## Version-dependent inference

The preprint used species-level Poisson GLMs and reported strong surface preferences for O. discolor and O. ferruginea, and a community exposure×surface interaction (p=0.003). The final paper instead used negative-binomial GLMMs with site intercepts and model selection. Its results state:

- Erythrodiplax funerea and Orthemis discolor: null models had the greatest weights (0.86 and 0.69).
- O. ferruginea: category model weight 0.36; the reported moderately high category contrast was p=0.08.
- Telebasis salva: surface model weight 0.38, green-surface contrast p=0.057, and the null model was not distinguishable from the surface model.
- Community models did not clearly improve over the null (several alternatives within delta AICc<2).

These are different analyses of the same programme, not independent confirmations. Do not carry preprint significance into the final-paper evidence table. Absence of a clear preference is not equivalence or proof that polarized surfaces have no effects. The final Methods list also contains duplicated/misnumbered model labels; the Results and model descriptions are retained without repairing the source.

## Consequence for the manuscript

This report informs the **boundary of the water-reflection hypothesis**: polarized surfaces can be experimentally studied, but this study does not test whether nighttime water reflection increases predators, reduces moth pollination or lowers wetland plant populations. It contributes zero ALAN lnRR effects and zero independent ALAN experiments. Keep it outside the three stage-specific pooled estimates. A second automated reviewer checked scope and version-dependent inference; this is not human double extraction or a raw-data reanalysis.
