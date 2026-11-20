# Burgazzi2026: aquatic movement and benthic standing stock

Primary DOI: https://doi.org/10.1002/ecog.07758
Data: https://doi.org/10.5061/dryad.rfj6q57k5

Methods, Results and author analysis code reviewed on 2026-11-20. See
[independent design review](../search/burgazzi2026_design_review.md).
This supplementary discovery is outside the existing indexed search export.
It does not increase that export's 124/190 reviewed count.

`source_manifest.json` records the original workbook's publisher-supplied hash,
not a successfully verified local workbook. Dryad metadata were accessible;
the original-data download returned API 401/public 403. No authentication
bypass was attempted. The author code was retrieved and checked against the
pinned commit URL in the manifest. External cache:
`work/gbif_annual/alan_experiment_sources/burgazzi2026/`.

No numerical effect was extracted. In particular, Table 2 coefficients for
square-root abundance are not log response ratios. Reported nonsignificance
was not encoded as a zero effect. The current main catalogue is unchanged.

Before numerical extraction:

1. Recover the original workbook and verify its recorded SHA-256.
2. Verify section, flume, header-block, date and subsample identifiers.
3. Resolve shared campaign/measurement status with Juvigny-Khenafou2024
   (10.1016/j.scitotenv.2023.168836); do not assume two independent programmes.
4. Fix the untreated comparator rule, retaining flume and header dependence.
   Separate final benthos from terminal drift and repeated-time supplements.
5. Keep derived downstream export separate from independent measurements.

These are movement and standing-stock responses. Predation, emergence,
pollination and plant reproduction were not measured. They cannot establish
the proposed predator-mediated wet-plant decline pathway.
