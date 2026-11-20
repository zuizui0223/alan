# Miller2017: non-warmed field stratum

Primary paper: https://doi.org/10.1098/rspb.2017.1195 . Original data: https://doi.org/10.5061/dryad.h7k26 ; recovered mirror: https://zenodo.org/records/4962533 . File MD5 `85935aca65b83672b6b48fd3266d2fed` matches Dryad. The copyrighted workbook remains outside the repository. `field_rows.csv` preserves its 48 field observations and adds a reconstructed block identifier.

## Allocation and scope

The [independent design inspection](../search/warming_predation2017_design_check.md) identified a missing warmed CM cell and duplicated warmed N cell in block4. These are different cages, not duplicated measurements. We retain both and do not guess corrected predator labels. `allocation_audit.csv` includes the empty factorial cell. The complete non-warmed stratum has24 cages, four blocks and three predator treatments. This subset supports three paired lighting effects and two dependent differences from the no-predator lighting effect. It does not estimate warming interactions.

Original codes are preserved. C7, CM and N are interpreted from the paper's treatment identities as Coccinella septempunctata, Coleomegilla maculata and no predator, respectively; an explicit workbook codebook was not found. Block=ceil(Row/2) reconstructs neighbouring row pairs and is not an original column.

## Calculation

Run `python database/alan_experiments/warming_sources/extract_field.py` from the repository root. Optional `--source` selects the external cache directory.

For each predator treatment, lnRR=log(mean lit/mean dark). The full six-arm covariance of block means is transformed with the log-ratio Jacobian. Predator contrasts are lnRR(predator)−lnRR(no predator); their two-by-two covariance preserves the common comparator and block pairing. Approximate intervals use t with3 degrees of freedom and the delta-method variance; four blocks provide limited precision. A scalar paired formula and block influence calculation verify the matrix calculations.

The C7 contrast is −0.90457, corresponding to a ratio of response ratios0.4047 (approximate95% interval0.0781–2.0966). The CM contrast is −0.10989, ratio0.8959 (0.1842–4.3586). Neither excludes1. These are final aphid population counts, not consumption or attack rates. They neither establish predator mediation nor demonstrate absence of an interaction.

`unheated_light_effects.csv` and `unheated_predator_contrasts.csv` are supplementary extraction products. Their covariance files must accompany quantitative reuse. Main catalogue inclusion and pooled estimates remain pending; these rows do not change the108-effect catalogue. The full-paper factorial test discrepancy is unresolved and was not used to reconstruct variance.
