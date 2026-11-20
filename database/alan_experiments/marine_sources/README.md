# Davies2015 marine board-level extraction: exposure hold

[Primary article](https://doi.org/10.1098/rsbl.2015.0080); [original supplementary data archive](https://www.ebi.ac.uk/biostudies/studies/S-EPMC4424621). BioStudies delivered all three supplementary files with byte counts matching its metadata. URLs and SHA256 hashes are retained. This resolves the earlier HTML response from the PMC attachment route. Source review: [Methods and Results](../search/marine_2015_review.md).

## Original-data calculations

The 1,692 records contain two panel observations for each of47 taxa on18 boards. No values are missing. Ten taxa use percentage cover;37 use individual counts. Cover is averaged and counts summed within boards, producing846 board/taxon measurements. Each control/low/high arm has six independent boards. Numeric suffixes are not treated as matched blocks.

All94 taxon-by-light contrasts are retained, including all boards and uncommon taxa.35 have positive arm means and finite lnRR;59 have at least one zero mean. No pseudocounts, deleted outlier or significance selection are used. Shared controls and cross-taxon covariance are retained;35 finite contrasts have covariance rank15 and are not35 independent studies.

Unstandardized mean differences are additionally stored for all94 contrasts, with joint covariance. Their units remain percentage points or individuals per board; these scales must not be pooled as one response. A comparison with zero sampling variance is marked in the audit, not given infinite weight. Differences are not presented as a substitute pooled ecological effect.

The paper's Methods describe dawn-on/dusk-off switching, conflicting with its nocturnal treatment description. Neither supplementary caption resolves timing. Until adjudicated, these estimates remain outside the80-effect catalogue. This is a marine colonization dataset, not a direct plant–pollinator or predator-rate experiment.

Reproduce with `python database/alan_experiments/marine_sources/extract_boards.py`. Original files remain in the external cache.
