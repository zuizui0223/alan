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
source completeness and the survival supplement require reconciliation.
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
parser. SupplementS1 cells remain unverified in this check. This is not a
completed human duplicate extraction or an explanation of the missing cell.
