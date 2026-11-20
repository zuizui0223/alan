# Retained-flower damage: paired uncertainty check

Primary DOI10.1016/j.isci.2025.112770; raw data https://doi.org/10.17632/z6xd3dnf7n.1 . This is an exploratory same-date paired subset, not an author-model reproduction or final meta-analysis adoption.

Independent AI recalculation from original `C_240124_TRAITS.csv` verified the raw SHA256 against repository metadata and reconstructed observed-flower means within plants, plant means within sites, then equally weighted paired site means. Author code lines113–114 confirm light labels; lines3530–3545 identify the stable pairs used here. Original data are retained externally; hashes and URLs are in `source_manifest.csv`.

| Species | Last date with both paired sites | Site pairs | RR | Paired delta lnRR variance | Approximate t1 95% interval |
|---|---|---|---:|---:|---|
| Leucanthemum vulgare | 2023-05-18 | MUSTELID, RAPTOR | 0.204393 | 0.391152173 | 0.0000723–577.681 |
| Malva moschata | 2023-07-13 | HYMENO, MUSTELID | 0.741273 | 0.000834282 | 0.513558–1.069958 |

`summarize_florivory.py` adds small-pair uncertainty to the existing extraction and preserves its hold statuses. The independent raw recalculation reproduced both ratios and variances. This is not human double extraction. No effects were added to pooled estimates or the108-row catalogue. Two species share MUSTELID at different dates, and missing cross-date covariance is not assumed zero.

Both intervals include1. With only two pairs, one pair's unusual value can dominate; the small variance for Malva does not establish population-wide precision. Centaurea jacea has one matched pair and is not assigned a replicated sampling variance. The latest supported date is not necessarily the last study-wide date for a species. No all-year temporal response was fitted.

Unobserved flower values are missing, not zero. For Leucanthemum at BADGER, one of25 plants has all flower values missing, so the flower endpoint averages24 plants. Undamaged observed flowers remain zero. Plant identities duplicated across sites in the wider source remain documented in `audit.json`; no silent global deduplication was applied. These endpoints concern retained sampled flowers, not lost flowers, whole-plant damage, observed slug-specific feeding, pollination or seed production. Existing source/support adoption holds remain.
