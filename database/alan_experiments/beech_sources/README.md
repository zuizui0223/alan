# Beech2025 original-workbook scope audit

[Primary-source review](../search/beech2025_review.md); [original data](https://doi.org/10.16904/envidat.628).

The original workbook has56 chemistry rows but855 leaf-damage rows from57 tree IDs. Its stated four treatment labels cover56 trees and840 leaves. The additional Alpthal63 tree has15 leaves under2200_100, a different spectral setting. The extra tree has no chemistry match. All rows are retained in the scope audit; none are merged into the4000_100 group.

`tree_support_and_spectral_scope.csv` links tree and leaf coverage; `installation_support.csv` counts trees/leaves by installation. The stated experiment has eight installation IDs, four per site, with seven trees per installation. This is not56 independent illumination replicates. Settings30/50/100 must not be relabelled as measured lux without the paper's calibration.

Leaf damage incidence and severity conditional on damage answer different questions. The initial audit computes support counts and descriptive tree damage only. The supplementary extraction below retains installation-level contrasts separately from the main catalogue. Potential Swiss programme overlap remains unresolved.

Run `python database/alan_experiments/beech_sources/audit_workbook.py`. The pinned original workbook hash is verified before reading.

## Supplementary installation-level extraction

Run `python database/alan_experiments/beech_sources/extract_installation_damage.py`.
The stated 840 leaves are aggregated into eight installations, one per dose per
site. The extra spectral treatment remains in the original scope audit.
Paper-calibrated nominal doses are 21, 35 and 70 lux; original `4000_*` labels
are retained (the paper describes the neutral-white lamps as 3700 K).

Three endpoints are retained: damaged-leaf proportion; mean leaf damage percent
including zeros; and mean damage percent among damaged leaves within each
installation. The second is not the ratio of total consumed area to total leaf
area. The third is not an equal-weight mean of tree-specific conditional means.
Each endpoint's installation values are averaged equally across the two sites.

`installation_effects.csv` contains nine dependent log response ratios and
paired-site delta variances. Approximate intervals use t with one degree of
freedom. `installation_covariance.csv` preserves all 81 covariance cells,
including shared dark installations and responses from the same leaves.
The covariance has rank one; nine contrasts are not nine independent studies.

At 70 lux the all-leaf mean ratio is 0.618 (approximate interval 0.411–0.930),
but the unusually small paired variance reflects agreement of only two site
scores. It does not establish robust precision or a general 70-lux threshold.
All other contrast intervals include one. The other all-leaf ratios are 0.442
at 21 lux and 1.149 at 35 lux, so the descriptive pattern is not monotonic.
No multiple-comparison correction or original GLMM reproduction is claimed.
No chemical mediation, predator mechanism or pollination effect is inferred.

[Independent raw-data check](../search/beech2025_extraction_check.md).
