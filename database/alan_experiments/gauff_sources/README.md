# Gauff2026: field seagrass–epibiota manipulation

Primary DOI: https://doi.org/10.1016/j.marenvres.2026.108273
Recovered institutional sources are linked and hashed in `source_manifest.json`.
The PDF is marked Journal Pre-proof; the inaccessible publisher version has
not been assumed identical. [Independent review](../search/gauff2026_epibiota_review.md).

Six field plots (three lit, three dark) each contain six quadrats. Three quadrats
per plot were cleaned of epibiota every two weeks. This is repeated field
removal, not only an immediate laboratory manipulation. Plot lighting alternates
spatially; random allocation is not documented. Quadrats, shoots and repeated
dates are not additional light replicates. The coastal marine setting does not
represent a replicated freshwater-wetland versus dry-land contrast.

Night-time field oxygen/PAM responses and subsequent standardized laboratory
daylight assays are different endpoints. Removal modified the reported
daylight-assay response, while night-time responses and biomass were not
detectably affected. The comparison does not measure predation or pollination.

## Source discrepancies retained

- Primary text assigns leaf-length p≈0.003 to lit versus dark cleaned plants.
  Supplement Table 11 instead reports that contrast as -24.084, SE13.311,
  adjusted p0.2574. Its p0.0027 belongs to uncleaned-lit minus cleaned-lit
  (+31.981, SE9.233). An illumination-by-epibiota model interaction remains
  reported in Table 10, but it is not the claimed pairwise light comparison.
- Oxygen Tables 2 and 3 have identical numerical/model entries below their
  differing net/gross-production headers. Do not choose one as correct.
- Table 7's header and Table 8's caption identify dark adaptation, whereas
  the primary Results assign these numbers to the simulated daylight assay.
- Supplement Table 1 is not represented as an ordinary populated text table;
  no numerical respiration result was inferred from its empty table object.

`reported_pairwise_comparisons.csv` preserves all 12 rows from Tables 8 and 11,
including SEs and adjusted p-values, without relabelling the ambiguous stage.
These are author-reported contrasts, not self-extracted arm lnRRs. Do not derive
lnRR from signed oxygen fluxes or sampling SE from rounded t values. Raw
plot/quadrat data, covariance and final-version resolution remain necessary.
The primary data-availability sentence contains no usable dataset identifier;
no public raw data were recovered in this pass. No author contact was made.

Run `python database/alan_experiments/gauff_sources/extract_reported_tables.py`.
External cache: `work/gbif_annual/alan_experiment_sources/gauff2026/`.
No main-catalogue or pooled effect was added.
