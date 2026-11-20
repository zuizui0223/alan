# Slug herbivory experiment: original field-data audit

2026-11-20. Sources: [primary article](https://doi.org/10.1016/j.isci.2025.112770), [author data v1](https://doi.org/10.17632/z6xd3dnf7n.1), [author analysis](https://github.com/GrognuzV/Raw-data-for-ALAN-x-herbivory). The Mendeley public API returned the original files and SHA-256 hashes; four downloaded field/site files passed those hashes. Laboratory growth, mortality and activity experiments were excluded entirely.

## Design and interpretation

The field experiment used illuminated/dark wildflower strips in Switzerland (six pairs in 2022, eight in 2023). Field herbivory endpoints examined here are from 2023. Potted *Centaurea jacea* were exposed for a day–night cycle; only damage accompanied by gastropod slime was counted. Wild *C. jacea*, *Leucanthemum vulgare* and *Malva moschata* damage was not exclusively attributable to slugs. The paper reports reduced nocturnal potted-plant damage and reduced wild-plant florivory severity. Predator avoidance is a proposed explanation, not a measured predator-treatment mediation effect.

This is useful counterevidence to a universal pathway from illumination to increased herbivory and plant decline: illumination can reduce a humidity-associated herbivore's damage. It does not test wet-versus-dry effect modification, explain the macro wet-affinity enrichment by itself, or measure reproductive benefit from reduced damage.

## What was extracted

`database/alan_experiments/slug_sources/extract_original.py` reads original files directly, preserves original missing observations, and exports source manifests, plant/period rows, site-date means, selected paired contrasts and same-support effect covariance. No missing site/species record is treated as zero. Means/SDs and covariance use independent sites/pairs, never leaves or plants as independent treatment replicates.

For wild plants, calculated endpoints are: proportion of observed organs damaged; mean damage across all observed organs; and conditional mean severity among damaged organs. Leaves and flowers are separate. The unconditional mean is an explicitly new original-data summary, not the paper's conditional severity estimand.

For pots, the author's 0.1 damage threshold and truncation of negative changes to zero were reproduced for the proportion increase; raw changes are retained alongside. Conditional severity changes follow the author's rule, with undefined values left missing. There are 49 negative proportion changes before truncation. A 24-hour aggregate was **not** constructed because interval linkage and plant identities require clarification.

## Source barriers

- Potted source: 416 period rows, 151 unique plant labels, 270 plant-label/date combinations, versus 216 exposed pots reported in the article. Required leaf-damage fields contain 1045 missing cells. Plant labels recur across dates; no duplicate label/date/period keys exist. These are not 416 independent plants, and a missing period cannot be filled with zero.
- Wild source: 550 rows and 549 labels versus 544 plants reported. Label `23-08600` appears at both MARTEN and HORNET on the same date. No silent correction or exclusion was applied.
- The source code changes pair mappings for three original pairs in later sampling. We used only the five pair identities that remain stable in the code; this is a documented support restriction, not reconstruction of the complete published analysis.
- Wild plants occur unevenly among sites and dates. The matched-date extraction has only one pair for *C. jacea*, so no replicated paired effect is produced for that species. Two pairs each support the latest common-date contrasts for *L. vulgare* and *M. moschata*.
- Most wild plants were measured once. “Latest common date” selects sampled individuals at the latest observed date within each pair; it is not a longitudinal final measurement of the same plant and not a common end-of-season date.

## Descriptive final-observed contrasts, not approved meta effects

The following are ln(mean lit / mean dark), from the explicitly restricted source support. They must remain in a holdout table pending source-identity and sampling-support review.

| Contrast | Pairs | lnRR | Paired delta variance |
|---|---:|---:|---:|
| Wild *L. vulgare*, leaf proportion damaged | 2 | −0.90887 | 0.08080 |
| Wild *L. vulgare*, flower mean damage, all observed | 2 | −1.58771 | 0.39115 |
| Wild *M. moschata*, leaf proportion damaged | 2 | +0.38145 | 0.09383 |
| Wild *M. moschata*, flower mean damage, all observed | 2 | −0.29939 | 0.00083 |
| Potted *C. jacea*, night proportion increase, 30 August | 3 | +0.39412 | 0.58555 |

The final observed night subset therefore does **not** reproduce the article's pooled-season negative direction; its uncertainty and small selected support must be retained, not hidden. For this subset, the lit/dark mean proportional increase is 0.09563/0.06448, with unequal within-site pot counts (1–7), but equal pair weights. This is not a contradiction of the fitted all-season estimate because the dates, support and adjustment differ.

The potted daytime ratio is undefined where the group mean is zero, so no pseudocount was added. Night conditional-severity variance degenerates to zero on two pairs because one pair has zero in both arms; it is explicitly flagged and cannot be used for inverse-variance weighting. Covariances are exported only for endpoints with identical pair support; no covariance is fabricated across differing supports.

## Decision

Original data are obtained and auditable site-level extraction is complete, but **zero effects from this extraction are approved for the master meta-analysis**. Resolve plant-label/sample-count discrepancies and specify the intended field estimand before promotion. Prefer an explicitly defined all-season site summary if the purpose is season-wide response; retain these final-observed contrasts as a separate estimand. If the protocol requires a literal common final time, wild cross-sectional damage cannot be relabeled to satisfy it.

Master tables unchanged; no commit. This source strengthens the need to separate reduced herbivory, altered visitation, pollen delivery and plant population change rather than assigning them one ecological “harm” sign.
