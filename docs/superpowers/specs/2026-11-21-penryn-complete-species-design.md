# Penryn complete-species composition and trophic contrasts

Continue the authorized ALAN biological-interaction analysis from main
`db82976f834cc048e76e69e6f621ba11e3ef33e5`. The new question is whether the
2016 biomass response of Lotus pedunculatus differs from three co-occurring
species, and whether the experimentally assigned trophic treatments modify
that difference.

## Design and choices

Use the four species with positive biomass records in all 54 physical
mesocosms and matching counts in original Supplement S1: Agrostis tenuis,
Anthoxanthum odoratum, Hypochaeris radicata and Lotus pedunculatus. Retain
printed taxon labels. Exclude Centaurea nigra because its complete biomass
table conflicts with S1. Do not fill absent observations with zero.

The primary estimand is Lotus ln(mean biomass lit / mean biomass dark)
minus the equal-weight mean of the other three species' lnRRs. Estimate
this within each of three trophic treatments and two light colours.
Keep the full covariance across species sharing mesocosms and colours
sharing controls. Summarize each light colour across the three trophic
treatments with equal weights. Estimate T2 minus T1 and T3 minus T2
differences explicitly; an observed light response alone does not identify
an herbivore or predator mechanism.

Alternative considered: a 13-species moisture regression duplicates earlier
branch diagnostics and mixes unequal survival support and unresolved counts.
Alternative considered: a four-species moisture slope conflates wetness with
one species. Use the named-species contrast instead. A per-mesocosm log ratio
is a secondary estimand sensitivity, not the same arithmetic-mean lnRR.

Intervals use the empirical covariance and conservative approximate t5
critical values. Also report Bonferroni intervals for the two marginal
comparisons and separately four trophic modifiers. These are post hoc,
conditional estimates, not a preregistered confirmatory test or a new
independent experiment. No pooled meta estimate is added. T1 includes
repeated insecticide/molluscicide treatment; T2 versus T1 is the whole
trophic treatment package. T3 versus T2 is predator addition. Species
identity, legume status and wet affinity cannot be separated here.

## Deliverable

Reproducible analysis, full covariance and contrast tables, source/selection
receipt, a readable figure and Japanese results note. Link the new analysis
from the current synthesis. Preserve existing 108-effect catalogue, three
meta-analysis outputs and macro estimates. Use a new draft PR; no merge.
