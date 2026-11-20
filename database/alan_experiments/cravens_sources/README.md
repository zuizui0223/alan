# Cravens2018: bat diet composition audit

Primary: https://doi.org/10.1111/1365-2664.13036
Original workbook and representative sequences: https://zenodo.org/records/4994918 (Dryad gm6kk mirror). Both files verified against publisher MD5; hashes retained.

Methods/Results read2026-11-20. Twenty locations were surveyed with nearby dark captures on first two nights and lit captures on third nights; this is not a simultaneous paired comparison. Six bat species are dependent outcomes of one programme. DNA occurrence does not measure prey abundance or predation rate. Reported453captures differs from297+151; Table3's188analysed samples matches the workbook exactly.

`audit_occurrence.py` implements the written sample/order occurrence definition, removing repeated detections of the same order within a sample. It separately reports sample prevalence and occurrence normalized across orders. Neither uses read depth as prey abundance. Under this sample-order definition MYLU declines. Follow-up reconstruction below resolves the direction through a different detection unit; preserve this initial output as a separate estimand.

Sample prefixes are not verified physical site IDs. No site-pair covariance or independent lighting variance was inferred. No main effect added. Independent primary-source review confirmed the experimental timing, measurement limits and188sample count.

## Resolved aggregation distinction

Follow-up2026-11-20 counted each detected prey taxon once within each sample, then summed detections across samples. This reproduces the reported grey-bat64.9% increase and positive MYLU direction. MOTU-key and full-taxonomy-key calculations agree, independently checked using an incidence matrix. Neither count is numbers of prey consumed. `moth_prey_taxon_descriptive.csv` reports all six species without a pooled estimate or inferred sampling variance.

The raw workbook contains490unique MOTU/taxonomy labels versus487reported in the article. Representative sequences number3078, matching the intermediate sequence-cluster count. No taxa were removed to force agreement. Therefore the grey-bat numerical agreement is a specific validation, not proof of full Figure1 reproduction. Site/night dependence and the remaining version discrepancy still prevent final meta inclusion.

Independent raw-workbook audit confirmed all six numerator/denominator pairs and ratios; retained490taxa and no exclusions.
