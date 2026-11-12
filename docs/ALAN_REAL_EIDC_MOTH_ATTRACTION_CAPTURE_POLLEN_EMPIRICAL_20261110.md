# ALAN — authentic light attraction versus ground-level moth capture and body-pollen
2026-11-10; original data / sampling schema audited in separate signed GitHub Actions.

## Sources and source fidelity

Original field data (not web-text values) were downloaded from the published [UK EIDC 2016 archive](https://doi.org/10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708). Four actual source CSV files were byte-checked against source-audited SHA256:
- F1_sites: 22 original rows, 20 valid matched source sites
- F2_moths: 1,599 individual captured moth rows, 20 sites, 609 transect and 990 light-trap individuals
- F3_overhead: 90 observations, 45 complete matched site-by-date lit/unlit survey sessions at 20 sites
- F4_pollen: 361 positive body-pollen record rows, 359 distinct source moth-reference strings; **one** reference not in original F2, two extra rows corresponding to duplicate-source reference, leaving **358 exactly linked moth IDs** (keep unmatched unknown)

[Actual authenticated results and stats](../results/alan_eidc_2016_real_moth_flight_capture_pollen_results_20261110_verified.json); [successful original source inference Actions](https://github.com/zuizui0223/alan/actions/runs/42767682403), artifact 13740117359. [Original F1-F4 relational keys audit Actions](https://github.com/zuizui0223/alan/actions/runs/42767214170).

## True matched light/dark source results

Each original pair represents the lit and unlit portions of one field-margin site, sampled at matching times. Multiple survey nights are NOT independent sites. Source matched pairs were retained without imputing an unrecorded moth capture event as an actual zero.

| Original outcome | Source total lit | Source total unlit | Equal-**site** mean log ratio exp-transformed | Nominal site-based t 95% CI | Source method limit |
| --- | ---: | ---: | ---: | ---: | --- |
| Overhead flight activity, sum of three passes (45 complete pairs, 20 sites) | 286 | 148 | **2.10554** | **[1.5970, 2.7760]** | 0.5-pseudocount log ratio; separate from paper's original mixed GLMM |
| Ground transect captured moths (45 observed-positive paired events, 20 sites) | 242 | 367 | **0.65798** | **[0.54135, 0.79974]** | F2 individual moth records do not establish a complete zero-capture effort ledger |
| Ground light-trap captured moths (40 observed-positive paired events, 20 sites) | 455 | 535 | 0.93721 | [0.77416, 1.13459] | The trap method itself attracts moths; method-dependent detection |

Nominal unadjusted two-sided site-mean t p values: overhead 1.955e-5, transect 2.507e-4, trap .4862. These three results are **post-hoc reanalysis of the SAME source publication, not independent evidence of three experiments or multiple-testing-corrected confirmatory tests**. The published mixed-model result overhead flight activity was about +70%; our crude raw pooled overhead ratio is 286/148=**1.932** whereas equal-site pseudo-count-adjusted ratio is **2.106**. Different estimands, neither supersedes the original paper.

Recorded F4 positive pollen source totals linked to F2 captured individuals: lit 137 genuine individuals, unlit 221 (of 697 vs 902 captured originals); recorded pollen-grain sums from matching references: lit **3,995**, unlit **13,386**. These totals are **not normalized for every sampling opportunity and are strongly count-weighted**. Do not call these successful stigma deposition, seed set or fruit production. Do NOT silently treat all F2 moths without F4 record as verified pollen negatives: that requires the supporting biological sampling protocol to establish which individuals were assayed.

An exploratory cross-site correlation between metered local lamp lux and equal-site flight ratio gave Spearman rho **−0.453**, nominal p **.0449**, n=20 matched sites. With only 20 sites, limited site ranges and many exploratory comparisons this is **not a robust dose-response effect** and must not be used to claim brighter streetlights repel moths.

## Interpretation

**Both apparently opposing effects occur together in the same original matched field dataset**: more flights at lamp height and fewer moths recorded on transects at ground level; a separate moth trap method gives an unresolved weak negative contrast. Positive moth aggregation cannot be equated with effective pollination or persistent floral-composition benefits. The result directly demonstrates why ALAN source-pollination claims must separate measurement height and trap mechanisms.

The 2016 survey is mainland South-East England, not the 970-island 2025 original PCNL pilot. It does **not establish** whether local streetlights increase moth flower visits on islands, nor whether 5-km dark-corridor connectivity changes that process. It is a valid independent moth outcome *benchmark*, not a same-location biological response for the island macro study.

## Next biological link gate

Require an actual nocturnal insect census or moth pollen-deposition/seed set at the **same island, site, time** as satellite exposure, or a defensible matched geographic independent census. The current ALAN source plant GloBI positive-visit fractions do not measure active moths. Do not merge the 2016 England field inference with 2025 island observational photo composition into one causal slope.
