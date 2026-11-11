# GloBI v0.11 × GloPL pollen-limitation experiments × historical nightlight

**2026-11-11 original-source empirical audit.** This is a direct exact botanical-name rejoin, NOT an ecological claim that moth pollination lowers pollen limitation, nor a causal test of ALAN disrupting effective nocturnal pollen delivery.

## Source contracts

- Complete interpreted GloBI v0.11, **24,453,954** original relationship rows (source SHA256 \`670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3\`). Previous [full-reproduction run #42924993614](https://github.com/zuizui0223/alan/actions/runs/42924993614) generated source-positive botanical name flags, saved under artifact 13815335900. 1,651 botanical names had positive moth-family flower-use source claims, 4,065 had butterfly-family claims. Positive historical global record, NOT current local effective pollination or no-claim biological absence.
- GloPL, Bennett et al. (2018/2019) **2,969 original pollen supplementation experiments / 1,265 published plant species** [Dryad DOI 10.5061/dryad.dt437](https://doi.org/10.5061/dryad.dt437). To avoid future light and false publication-year matching, this **historic ALAN subset** is the **1,419 actual original experiments conducted in 1993–2015** with verified actual experiment year, original location, independent publication study key and complete pre-experiment 1992–2014 Li–Cao PCNL. It covers **341 original study keys, 507 original accepted plant names and 597 source geolocations**. Signed GloPL+PCNL joined member SHA256 \`a2b65a9a90e1cb97ed6ee366da8e8b86e67d14b51ded918284eaec3d917174a3\`, [original run #42754948691](https://github.com/zuizui0223/alan/actions/runs/42754948691), artifact 13735831120. The original 2,969 rows were not all given a historical light estimate, and this 1,419 subset cannot describe all GloPL experimental populations.
- **Exact scientific botanical accepted names** joined by equality; 36 original GloPL experimental records had 15 infrataxon/unresolved names and were retained as unmatched source-record-unknown (not assumed true non-moth); 2026 GloBI names do not establish recorded field interactions in 1993–2015.
- Measured PL response is the published **log response ratio of experimental supplemental pollen reproductive output versus the naturally pollinated control**, not a general flower abundance, field moth census, or nocturnal-specific pollinator-deficit estimate.

## True exact-source botanical overlap

| GloBI documentation state of original GloPL experimental species | GloPL experiment records | Distinct plant names | Independent GloPL study keys |
| --- | ---: | ---: | ---: |
| Moth-positive full source | **170** | **41** | **42** |
| Moth-positive excluding iNaturalist interaction sources | 100 | 26 | 26 |
| Moth-positive excluding iNaturalist and trophiCH source aggregators | 97 | 23 | 25 |
| Butterfly-positive full source | 248 | 69 | 59 |
| **Both** moth and butterfly positive documented | 144 | 31 | 32 |
| Moth recorded, butterfly **not recorded** | 26 | 10 | 12 |

Moth-only **recorded** is NOT proof that butterflies never pollinate. Among 507 original source plant names, 41 with documentary positive moth visits is a lower bound of documented source coverage, NOT a denominator for biological nocturnal-pollination prevalence.

## New post-hoc empirical pollen-supplementation comparisons

Fit original source **study-equal weighted WLS** with publication/study-cluster t CI, controlling year of actual experiment, verified mean historical pre-experiment log1p(PCNL), latitude and longitude with quadratic terms and the original experimental measurement/augmentation categories.

| Outcome and contrast | Beta in log reproductive-output ratio units | Nominal 95% study-cluster CI | Nominal p |
| --- | ---: | --- | ---: |
| Pollen-limitation experimental effect, GloBI moth-positive source taxon versus **no recorded moth claim** | **−0.2104** | [−0.4120,−0.0089] | .0408 |
| Same association using non-iNaturalist GloBI moth-positive records | −0.2470 | [−0.4845,−0.0095] | .0415 |
| Further exclude trophiCH source | −0.2283 | [−0.4725,+0.0159] | .0668 |
| Moth-positive × actual historical pre-experiment ALAN light | **−0.0344** | [−0.1595,+0.0907] | .589 |
| Moth-positive source effect adjusting for butterfly-positive recorded status | −0.2340 | [−0.5768,+0.1089] | .180 |
| Only **12 genera** represented by BOTH moth-documentation states (141 original experiments, 34 GloPL studies, genus FE) | **+0.1212** | [−0.4021,+0.6446] | .641 |
| Only **14 studies** comparing BOTH documented and undocumented species in the same study (278 experiments) | −0.2766 | [−0.9104,+0.3572] | .363 |
| Independent previously fitted complete 1,419 case GloPL outcome × historical ALAN, without source-positive GloBI term | +0.0107 | [−0.0468,+0.0683] | .715 |

For five overlapping source-label main contrasts, **Benjamini–Hochberg q for the all-moth group = .0692**, not q less than .05. These outcomes are **exploratory after multiple specification checks**, not five independent predictions or a discovery. The apparent weaker PL in moth-record-positive plants is not robust to source family correction, genus-within comparisons, same-study comparisons or requiring an independent butterfly-positive documentation control.

**Most defensible result:** all 2,969 published GloPL cases measure pollen-supplementation response, but only 1,419 have this exact historical source exposure. Within that subset, 41 botanical names have documented GloBI moth links. A suggestive *cross-study documentation-group* PL difference exists, **but GloBI moth-source status does not show a differential ALAN-related pollen limitation response**. This supports data overlap and a conditional research question, NOT evidence for or against local loss of effective moth pollination.

## Reproduction and where to go next

A self-contained source-locked Python reproduction was supplied in the current-turn user-facing code/result ZIP; it reads the two exact GitHub Actions ZIPs, independently verifies original source SHAs, checks input row and distinct-name identity, and recomputes WLS/study-cluster CI, source-exclusion flags, FDR and genus/study sensitivity. **Only public aggregate result JSON** is committed here: \`results/alan_GloBI_v011_GloPL_PL_priorALAN_real_crossjoin_20261111.json\`; no original species×experiment coordinates, private plant-photo identity ledger or original third-party rows are published.

Two independent next steps are scientifically justified: (a) expand the source-only GloBI overlap to the *full* 2,969 GloPL experiments WITHOUT pretending all of them have historical PCNL; (b) identify actual adult nocturnal moth effect on seed set via independent day/night exclusion and wild-flora local-study evidence, then test historical ALAN only in populations with matched source effort and verified observation year. No existing fd/island Chapter 1 result or other ALAN PR is overwritten. PR #40 remains a Draft.
