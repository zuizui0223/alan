# Cravens2018: bat diet composition audit

Primary: https://doi.org/10.1111/1365-2664.13036
Original workbook and representative sequences: https://zenodo.org/records/4994918 (Dryad gm6kk mirror). Both files verified against publisher MD5; hashes retained.

Methods/Results read2026-11-20. Twenty locations were surveyed with nearby dark captures on first two nights and lit captures on third nights; this is not a simultaneous paired comparison. Six bat species are dependent outcomes of one programme. DNA occurrence does not measure prey abundance or predation rate. Reported453captures differs from297+151; Table3's188analysed samples matches the workbook exactly.

`audit_occurrence.py` implements the written sample/order occurrence definition, removing repeated detections of the same order within a sample. It separately reports sample prevalence and occurrence normalized across orders. Neither uses read depth as prey abundance. Under this definition MYLU declines, whereas the paper describes an increase; therefore these outputs are a denominator audit, NOT reproduction of Figure1 or validated meta effects. Resolve taxonomic collapsing and figure aggregation against original code before extracting a contrast.

Sample prefixes are not verified physical site IDs. No site-pair covariance or independent lighting variance was inferred. No main effect added. Independent primary-source review confirmed the experimental timing, measurement limits and188sample count.
