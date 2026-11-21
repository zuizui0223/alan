# Penryn species biomass: source rows and conditional descriptive contrasts

2026-11-21. Original author data: https://doi.org/10.6084/m9.figshare.20304039,
PDF https://ndownloader.figshare.com/files/36265164. Article:
https://doi.org/10.1002/ece3.9157 (PMC9352868). The cached primary XML
sections 2.2 and 3.1 were read for tissue and survival scope; the public PMC
route subsequently returned a challenge, not an additional successful read.

The final August2016 harvest dismantled the mesocosms and separated roots
and shoots. Species total biomass includes both organs. This differs from
the clipped aboveground material in the Lizard grassland programme; these
two catalogue endpoints must not be silently pooled as identical measures.

The article analyses 13 species after excluding five rare species, while
community biomass includes all surviving species (up to17 per mesocosm).
The printed data yield665 rows for13species over54mesocosms. Species-table
sums do not reproduce all community totals: the largest deficit is64.20g
at A6. Missing rare species are a possible explanation, not identified
values to be redistributed. The original community table is unchanged.

Five species have records in all54mesocosms: Agrostis tenuis,
Anthoxanthum odoratum, Centaurea nigra, Hypochaeris radicata and
Lotus pedunculatus. Names are kept as printed, without an unrecorded
taxonomic substitution. Other species have fewer rows. No unlisted
species-by-mesocosm entry is filled with zero, nor interpreted as mortality.

There are116 observed species-by-light-by-trophic cells, each with3–6
recorded mesocosms. Prunella vulgaris L2/T1 has no parsed source rows.
This conflicts with a blanket reading of the article's claim that all13
analysed species have at least three surviving mesocosms in every cell;
the subsequent TableS1 check below resolves the exclusion annotation, not the missing biomass value.
Of78 planned lit/dark comparisons,77 yield finite descriptive point ratios;
the remaining comparison is explicitly held, not assigned ratio zero or one.

Festuca ovina E2 has total0.20g, root0.15g, shoot0.50g. Visual inspection
of original PDF page6 confirms the discrepancy is in the printed table.
No decimal correction is invented. The recorded total is retained and its
L3/T2 descriptive comparison flagged; the organ sum is not substituted.

The new outputs preserve source numeric strings, page references, treatment
allocation checks, arm support and every planned comparison. Means are over
recorded species rows. No confidence intervals, cross-species covariance,
wetness regression, trophic mediation or pooled effect is supplied.
The13species share experimental units and controls; none adds an independent
programme. The108-row main catalogue and114-row reported-statistics table
are unchanged. These contrasts do not reproduce the original transformed
mixed-model results or estimate unconditional population biomass/survival.

Run `python database/alan_experiments/plant_sources/extract_species_biomass.py SOURCE_DIRECTORY`.
`species_biomass_audit.json` records the PDF hash and checks. The manuscript
table is [the complete species-by-treatment display](../../../manuscript/biotic_synthesis/TABLE_penryn_species_biomass_ja.md).

An independent AI reader re-parsed the original PDF and reproduced all77
finite point ratios. All38 Prunella occurrences in the species table were
accounted for; L2/T1 is absent from that table text, rather than lost by the
parser. That first check did not verify SupplementS1. The subsequent visual reconciliation below supersedes that access limitation; neither check constitutes human duplicate extraction.

## SupplementS1 reconciliation (2026-11-21)

Visual transcription of original supplement page5 covers18species×9cells=162counts. TableS1 describes mesocosms with surviving individuals, not individual survivor counts. All printed row totals were checked. The162cells are preserved in `species_survival_reconciliation.csv`; the13species biomass comparison covers117cells.

Of these117cells,113agree directly. Prunella vulgaris A/T1 has1 surviving mesocosm, with an asterisk explicitly excluding that sample from its GLM; the biomass table has no corresponding row. Thus the omission is consistent with a documented author exclusion, not extinction. S1 total39 vs biomass38 is explained at the count level, but no biomass value or mesocosm ID for that excluded observation can be reconstructed.

Centaurea nigra remains inconsistent: S1 C/T1=5, A/T1=4 and W/T3=5, while the biomass table has6 positive rows in each cell. S1 total50 contrasts with54 biomass rows. The other11species match all nine cells. Consequently, five species have complete54-row biomass tables, but only four of those also have54 surviving mesocosms confirmed by S1. The Centaurea T1 amber/dark, T1 white/dark and T3 white/dark descriptive ratios depend on discrepant cells; they are not cleared for synthesis. No values are removed or replaced to force agreement.

Run `python database/alan_experiments/plant_sources/reconcile_species_survival.py SOURCE_DIRECTORY`. The supplement hash is saved. An independent AI reader checked the13species×9cells from the rendered table and confirmed the same4differences (one documented exclusion, three unresolved count conflicts). No survival effect, new independent programme, pooled effect or additional catalogue row is inferred. Festuca's organ/total discrepancy remains separate and unresolved.
