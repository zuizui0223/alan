# Amphibian2018: original data recovered; field phase reviewed

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6053935/) field Methods/Results read through indexed primary text; indoor Methods review remains incomplete. Initially40 outdoor mesocosms formed ten lighting groups. One illuminated group was removed, leaving four lit/five natural-light groups, each with four lid×toad treatments. A lamp served a group, not a single mesocosm. Lids altered both colonization and received illuminance, so their interaction is not a pure colonization mechanism.

Field metamorphic duration and relative periphyton fluorescence declined; survival and metamorphic mass showed no treatment response. Predator colonization is distinct from predation. Indoor juvenile effects must not be added as an independent field experiment.

[Dryad](https://doi.org/10.5061/dryad.pq54gh5) and [Zenodo](https://zenodo.org/records/4932781) identify the same workbook; bytes and MD5 agree. No raw effects extracted yet. Own audit in `amphibian_sources` identifies differing group labels across sheets and a long-format calendar-year discrepancy. Confirm physical lighting layout before estimating group-level covariance. This is partial phase-specific review, not a completed full-paper eligibility decision.

Supplementary DOCX was recovered from Figshare article6668306. It confirms shared lamps but provides no mesocosm-ID layout resolving sheet labels. Its title says juvenile growth increases, unlike the main article; retain this source-version discrepancy.

## 2026-11-20: group partition resolved from original workbook

The 14 label differences are not a partition contradiction. All 18 by 18 shared-mesocosm pair comparisons agree on group membership. The nine block names are a bijective relabelling: A:A, B:F, C:B, E:C, F:G, G:D, H:H, I:E, J:I (toad sheet : colonization sheet). `lengthofmeta` uses lowercase versions of the toad-sheet labels. Original labels are retained.

Mesocosm IDs join the water-quality and invertebrate sheets to the colonization sheet one-to-one for all36 mesocosms; light assignments agree. Metadata explicitly defines Mesocosm as an individual identifier and Block as a spatial group of four mesocosms. This supports using the source-labelled block mapping for those response tables. An independent reader also confirmed lid/toad treatment agreement for these36-ID joins. The audit exports `verified_block_label_map.csv` and `verified_mesocosm_join_support.csv` and checks the full324-element shared-membership matrix.

This removes the sheet-label obstruction; it does not independently verify the physical lamp layout, fix the juvenile calendar-year discrepancy, or finish indoor Methods review. No effect was extracted and primary-review counts are unchanged. The previous request for a layout to resolve mere label differences is superseded by this data-level check; a physical lighting-design check remains separate.
