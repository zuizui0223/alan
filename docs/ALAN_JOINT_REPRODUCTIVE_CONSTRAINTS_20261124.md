# Joint reproductive pathway constraints, wet ALAN cohort (2026-11-24)

Exploratory reanalysis of prior Pladias species-account extraction (wet_reproduction/species.csv, blob c603e16c96f5ee90d06df05205836dba9d206198) joined on sp to frozen b780370 wet_network/species.csv (blob 388a950643cda5769a41fbcecc8b4d4194448f7c), with exact accepted label and candidate status checked. Excluded all duplicated accepted-name rows and retained verified accounts: 270. Original source remains Pladias, Chytrý et al. 2021 DOI 10.23855/preslia.2021.001, generative reproduction classification Chrtek 2018; pollen vector and reproduction type from original account sources. No new primary experiments or ALAN effect sizes.

Joint proxy = insect pollination recorded without wind or water pollination (selfing/geitonogamy records allowed) AND exact allogamy or allogamy self-incompatibility label AND only by seed/spores. A positive proxy is not measured obligate insect dependence, lack of actual selfing, or lack of demographic persistence. Mixed mating codes not treated as strict allogamy; SI proxy is exact Pladias code, not experimentally verified genetic SI. Empty source values are missing. Seed-only label is species-level classification, not observed recruitment. Unknown in any component remains unknown, including logical AND. Record counts are not independent experimental studies.

| proxy | candidate yes/known | reference yes/known | strict matched difference (percentage points) | nominal permutation p |
|---|---:|---:|---:|---:|
| insect | 22/46 | 99/198 | -0.864 | 0.9345 |
| strict | 10/45 | 30/188 | -0.136 | 0.9864 |
| SI | 4/45 | 12/188 | -5.909 | 0.4205 |
| seed | 8/49 | 31/207 | 3.514 | 0.6311 |
| both | 4/44 | 18/181 | -2.083 | 0.7667 |
| joint | 1/44 | 4/181 | -1.504 | 0.7818 |
| siconstraint | 1/44 | 1/181 | 0.000 | 1 |

Matched strata = major site region, survey-site-number band, island-share band; require both candidate and reference. Equal candidate weighting of within-stratum differences; 9999 label permutations within strata, seeded LCG 7125. Unadjusted and region-only sensitivity in JSON. Nominal, exploratory, no multiplicity or phylogenetic/site-overlap correction. Joint full cases 44/181; strict matched 41/104 across 16 strata. One candidate (Ranunculus sceleratus), four reference labels (Frangula alnus, Salix triandra, S. purpurea, S. myrsinifolia) have joint proxy; this is not a difference in measured reproductive response. Four-profile counts candidate 1/3/7/33, reference 4/14/22/141. No clear enrichment. Cannot infer equivalence, no light-induced causal mechanism, wet main result unchanged.

## Broad-source acquisition audit

Ferrer et al. online 2024, issue 2025, DOI 10.1093/aob/mcae056, PMC11805948. 5,686 hermaphroditic taxa; 1,160 reported SI includes 725 inferred from congeners and 435 direct prezygotic evidence. Keep SS, PSC, SC and confirmed SI distinct. Original supplemental file https://pmc.ncbi.nlm.nih.gov/articles/instance/11805948/bin/mcae056_suppl_supplementary_data.xlsx served preparation page rather than binary; local environment DNS disabled; prior island repository workbook hash 8ef2f5ac99780ca19a15b847442272457634064577f8a12fbb9710f5521e5899; main References and SI Phenotype data sheets. Checked island curated evidence blob 2082c72531cdaede565df2f9edeb4bb84da8ac78: only three Ferrer synonym-recovery rows, none exact matching wet cohort. Not entire Ferrer source.

Fonseca et al. 2026 Dryad 10.5061/dryad.73n5tb3c9: German native angiosperm study (1,884 species overall, not 1,884 heterostyly measurements). CC0 data.zip direct file stream https://datadryad.org/downloads/file_stream/4700463 returned 403. README verifies Data_OtherSources_Heterostyly.csv semicolon-delimited, Species, Polimorphism (0 explicit homostylous / 1 heterostylous), Hercogamy, Reference, doi, Obs; several other traits require separate database access. Check explicit species-level evidence versus genus-level coding, conflicts and missingness; no unrecorded=0. Simon-Porcar et al. 2024 DOI 10.1038/s41467-024-45118-0 provides a 247-genera review of style-length polymorphism, not a species-level absence database. No broad SI/heterostyly coefficients reported.

Source acquisition unresolved; use original supplementary workbook and CC0 zip when available. Do not treat previous evidence from the few focal species as a denominator for all wet species. The present analysis uses existing accessible species traits, not those two new DBs. Calculation executed in-session; standalone Node replay not independently run in Actions. Main and main manuscript unchanged.
