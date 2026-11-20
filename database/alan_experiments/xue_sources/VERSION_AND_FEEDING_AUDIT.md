# Xue publication versions and feeding-data audit

The [2025 preprint](https://doi.org/10.64898/2025.12.01.691521) and [2026 article](https://doi.org/10.1002/ecy.70458) have the same ordered five authors, focal organisms, station, 2023 experiments, monoculture and community sample sizes, and the same 48.0%/87.2% community contrasts. Treat them as one publication family for extraction; this is a manual content adjudication, not an asserted DOI redirect. The bioRxiv API still reports `published=NA`. Separate report identifiers are retained. The final article remains the extraction source; no second copy of its four biomass effects is added.

Preprint Methods 2.2–2.3 confirm three light and three dark monoculture zones and zone/bucket nesting for feeding models. They do not resolve community bucket-to-light mapping. Reading the preprint therefore does not release the community hold.

## Original workbook audit

`audit_feeding.py` validates the previously downloaded workbook hash and retains assigned measurements. The 100 buckets yield 160 intended plant measurements: 157 observed, three missing. Another 240 blank cells represent plants not assigned to a no-choice bucket. Those structural blanks are not failed observations. Missing assigned observations occur at bucket48/Vallisneria, bucket99/M. spicatum and bucket66/M. aquaticum.

Among 157 observed initial-minus-residual mass differences, 83 are negative. Residual mass plus reported difference reproduces a constant initial estimate for each plant (0.453, 0.562, 0.430 and 1.212 g). These are signed estimates, not directly counted tissue consumption. Preserve negative values; do not set them to zero or compute a consumption lnRR. Zone support differs among assay/species/snail/arm cells. Choice-assay plants share buckets and illumination zones.

A mean-difference analysis would require explicit treatment of shared zones, initial-mass estimation and assay differences. The current audit supplies support tables, not an additional meta-analysis effect or evidence that negative values represent feeding below zero.
