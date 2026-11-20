# Lake2026: sequential illumination and zooplankton redistribution

## Primary review (Methods and Results read)

[Published PDF](https://www.nature.com/articles/s41598-026-60447-4.pdf), obtained2026-11-20: one anchored platform illuminated Lake Roś with a70W HPS lamp3m above water. Two2017 campaigns each repeated Night–Day–ALAN–Day three times. These are sequential observations at one installation, not simultaneous randomized illuminated/dark sites. Two hauls per depth were subsamples. PAR was measured in µmol photons m⁻² s⁻¹; approximate lux conversions are not independent measurements.

Fish echograms were opportunistic, qualitative illustrations, explicitly excluded from quantitative analyses. They do not supply fish-abundance or predation-rate effect sizes. Zooplankton outcomes distinguish mean vertical depth from integrated density within0–6m. Nocturnal deepening was supported for selected taxa and gravid females; consistent density reductions were not general. Chaoborus depth values included reconstructed sediment-associated positions, and the taxon was excluded from density models.

[Zenodo19249192](https://doi.org/10.5281/zenodo.19249192) supplies an Excel workbook. Hold main replicated-experiment extraction: temporal order/carryover, date linkage and shared infrastructure require resolution. A separate within-site temporal contrast could describe redistribution; it cannot establish consumptive loss or provide multiple independent lake replicates.

## Acquisition and independent archive inspection

- Cached published PDF:2,747,911bytes; published text extracted locally. Unlike the indexed HTML early-version landing page, the PDF identifies Scientific Reports16:30965 and contains complete Methods/Results.
- Cached preprint PDF separately; no replacement of the published Methods by preprint text.
- `ALAN data set.xlsx`:42,492bytes; MD5 `037689309905516e3f9a291bd8e6b560` matches Zenodo metadata. Access URLs/statuses are recorded in sibling cache `lake2026/access_audit.json`; metadata saved in `zenodo_metadata.json`.
- Workbook sheets: MVD192rows, Integrated-density192rows, Body-size432rows. MVD/density have eight taxon labels across May/Jun; each taxon/month contains DAY6, NIGHT3, ALAN3 rows. These counts match the differing day/night schedule, not24 independent illumination installations.
- Available keys are Month, Species, Treatment, n (and Depth in body-size). No explicit calendar date, installation ID, raw depth-specific haul count or time-series ordering table is supplied. Pairing NIGHT n1 with ALAN n1 cannot be asserted from the integer alone.
- MVD has97 missing gravid entries; density has96. Body-size has10 missing values. These are retained as missing, not zeros.
- The density sheet includes Chaoborus entries although the published density model excludes that taxon. Follow endpoint eligibility rather than importing every numeric cell into the main catalogue.
- Another reviewed Lake Roś report (PMC10997633; `fulltext_batch4.md`) also uses one boat light and repeated nights. Actual date/platform/data overlap with this publication remains unverified; do not automatically count it as an independent programme.

## Extraction boundary

No effect was computed in this check. Preserve PAR rather than inventing exact lux, do not digitize qualitative fish echoes as calibrated counts, and keep MVD in metres separate from density. Even a precise temporal standard error would not recover independent installation replication or remove fixed treatment-order confounding.
