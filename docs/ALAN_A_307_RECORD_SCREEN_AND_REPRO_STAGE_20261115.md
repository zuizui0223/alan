# ALAN A: full-paginated source discovery, title triage, and reproductive effect compatibility (2026-11-15)

## Actual search executed

GitHub Actions SUCCESS #43692485024: https://github.com/zuizui0223/alan/actions/runs/43692485024

OpenAlex title/abstract searches: 20 queries, all reported pages retrieved (each up to 20 pages at 100 per page), all queries succeeded, zero reported API failures. Total matched search rows 460; DOI-first deduplicated records 307 (previous first-page run 279). These 307 are bibliographic candidates, NOT independent experimental studies, fecundity effects, or a completed PRISMA census.

Source and triage artifact #14152429792 (ZIP SHA256 cd7cbcb8737732b26dec550fd5231a97d9e94a574ef666ee90185b8328d235f6): https://github.com/zuizui0223/alan/actions/runs/43692485024/artifacts/14152429792 . It contains source index CSV, the title-level triage CSV, and both audit receipts. All 307 titles retained.

Automatic title queue counts (NOT final inclusion/exclusion): registered reproductive papers 8; registered observational gradient 2; registered visits/phenology mechanism 2; other plant or pollinator title 52; review, thesis or data record 27; unresolved title 216. The 28 extra records gained by pagination yielded no newly certified primary reproduction experiment. Title screening cannot exclude or certify anything: no full-text decisions have been made in this automatic queue.

OpenAlex exact DOI recall among the 14-entry source registry: 12 retrieved, two missing despite known independent evidence. Missing Boom et al. 2020 DOI 10.1038/s41598-020-75471-1 (genuine manipulated floral reproductive study); missing Soteras et al. 2025 DOI 10.1016/j.ufug.2025.129070 (observational urbanization gradient). Hence complete pagination is not a complete literature census. Existing reproductive source registry: 12 publications (8 manipulated, one mixed-design, three observational), not 12 certified independent experiments.

## Reproductive values recovered and prevented from being falsely pooled

| Study | Verified form from original or author source | Source-based conversion | Why it is not presently a common A log-RR |
|---|---|---|---|
| Li et al. 2026 (10.64898/2026.07.10.737511) | 20/27 light array vs 8/26 dark PLANTS with at least one seed; 3 real experimental fields, 3 pooled species | (20/27)/(8/26) = 2.4074; log ratio +0.87855 | Plant seed occurrence, NOT fruiting among originally marked flowers; plant-level independent SE not appropriate with clustered treatment and species |
| Lu et al. 2026 JPE (10.1093/jpe/rtag160) | Pedicularis kansuensis seed formation 56% lower as reported by authors | light/dark factor 0.44; log factor -0.82098 | No recovered per-cell sample sizes, outcome denominator or uncertainty; Euphrasia regelii non-significant effect is not a true zero |
| Giavi et al. 2020 (10.1038/s41598-020-68667-y) | Epilobium angustifolium logistic model beta -1.37, SE 0.726 for selected capsules producing seeds, illuminated vs dark | odds ratio exp(-1.37)=0.2541, Wald 95% CI 0.0612-1.0544 | Odds ratio is NOT RR; estimated rate conditional on capsule sampling, meadow-level treatment clustering |
| Boom et al. 2020 (10.1038/s41598-020-75471-1) | Author-defined estimated seeds per originally exposed flower for green light vs dark | estimated ratio about 0.515, log ratio about -0.6636 | Model product of fruit set, uninfested fruits and ovule fertilization, not observed final seed total; overall estimated seeds ANOVA p=.63 |
| Macgregor et al. 2019 (10.1002/ecs2.2550) | Full-night reproductive-success direction positive, original lighting x day/night access interaction reported non-significant | No reconstructed source-level interaction beta or SE | Author R scripts accessible but 3 source CSVs unavailable; not a precision-zero finding |
| Knop et al. 2017 (10.1038/nature23288) | Published 13% fruit-set decrease | Keep as authors' reported wording, no unverified coefficient | Need actual denominator, original independent treatment units, exact SE before assigning a universal log-RR |

The effects above have different endpoints, scales, control designs and independence structures. A common global log-RR or night-dependence slope has NOT been estimated. Source-computation script scripts/analysis/alan_A_effect_type_gate_20261115.py and its six tests protect these distinctions.

## Critical implementation clarification: night-dependence is not necessarily from a different PAPER

A valid V can use separate daytime-only and open-access experimental observations at DARK treatment sites from the same paper as the light comparison, if measured independently from the ALAN response. In that case V = 1 - R(dark,day-only)/R(dark,open), and y = log(R(lit,open)/R(dark,open)). They SHARE dark-open controls, creating statistical covariance; source treatment cells and covariance are necessary for an unbiased slope fit. This clarification does not authorize using illuminated reproductive outcomes to invent V or using a different country's species-level V as though it were site/year-matched.

More fundamentally V is nighttime-access dependence, NOT automatically light-sensitive insect dependence: without proven nocturnal insect effectiveness and light response, label the inference accordingly. Nighttime vertebrates can provide effective pollen. Keep direct plant physiology, photoperiod, flowering phenology, light spectrum and light geometry as competing mechanisms rather than assuming a moth-only explanation.

## Gates and outcome

Completed: OpenAlex pagination and full preservation of 307 titles, title-based retrieval prioritization, source-endpoint compatibility rules.
Not completed: full-text independent screening of 307; matching same population/experimental dark-exclusion dependence to lit reproductive effects and their covariance; a meta-regression beta_V and CI; island-vs-mainland experimental interaction.

Next primary sources with highest inferential leverage: Macgregor 2019 unpublished/missing numerical original CSVs, Lu 2026 ALAN x four pollinator-access arms x two plant species source cells, Knop 2017 original plant fruit outcomes, Boom 2020 original light-allocation-level fruit and infestation data. All post-2026-11-15 searches must report their date and deviations: this is a protocol frozen AFTER seeing prior study outcomes, NOT retrospective preregistration.
