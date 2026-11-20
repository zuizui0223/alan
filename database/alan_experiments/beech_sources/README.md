# Beech2025 original-workbook scope audit

[Primary-source review](../search/beech2025_review.md); [original data](https://doi.org/10.16904/envidat.628).

The original workbook has56 chemistry rows but855 leaf-damage rows from57 tree IDs. Its stated four treatment labels cover56 trees and840 leaves. The additional Alpthal63 tree has15 leaves under2200_100, a different spectral setting. The extra tree has no chemistry match. All rows are retained in the scope audit; none are merged into the4000_100 group.

`tree_support_and_spectral_scope.csv` links tree and leaf coverage; `installation_support.csv` counts trees/leaves by installation. The stated experiment has eight installation IDs, four per site, with seven trees per installation. This is not56 independent illumination replicates. Settings30/50/100 must not be relabelled as measured lux without the paper's calibration.

Leaf damage incidence and severity conditional on damage answer different questions. This audit computes support counts and descriptive tree damage only; it adds no effect or sampling variance to the80-effect catalogue. Shared installations, source model definitions and potential Swiss programme overlap require adjudication before integration.

Run `python database/alan_experiments/beech_sources/audit_workbook.py`. The pinned original workbook hash is verified before reading.
