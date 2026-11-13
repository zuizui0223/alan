# ALAN — original four-arm day/night reproductive complementarity audit

**Authentic source replay SUCCESS: 2026-11-14 JST.** [Real original-source GitHub Actions #43456561380](https://github.com/zuizui0223/alan/actions/runs/43456561380). [Aggregate-only artifact #14054815172](https://github.com/zuizui0223/alan/actions/runs/43456561380/artifacts/14054815172). [Draft PR #58](https://github.com/zuizui0223/alan/pull/58).

## Scientific question and direct experimental outcome

Kendall & Nicholson, Ecology Letters (2025), DOI [10.1111/ele.70036](https://doi.org/10.1111/ele.70036), originally studied day-only, night-only and unrestricted open pollination, including complete visitor exclusion. The authors have **already published** day versus night and day/night versus open comparisons and discuss diel complementarity and possible functional redundancy. Our incremental original-data analysis quantifies **four-arm factorial nonadditivity** at true original field site×experiment year×season×botanical species×fruit-or-seed response. It does NOT newly discover the basic importance of open-pollinated controls or experimentally manipulate nightlight.

Source-verified original dataset: 1,652 original treatment rows from 140 papers, author Git blob SHA 39a82d966afcdde670603add1c881ace8791ffe9. Real original field-date/site day/night source: 301 paired fruit/seed outcome comparisons, 118 independent original studies, 123 original plant species. Source individual original study matching preserves all original published field site, actual experiment start/end year, monthly season, output metric and measurement units; day/night weighted treatment means for ALL 301 original paired cohorts were reproduced to original floating-point tolerance. Original source previous satellite Li–Cao PCNL exposure includes 30 MD5-audited annual GeoTIFF years 1992–2021, sampled at the true field site only strictly before the experimental year.

Of original 301 paired day/night reproductive outcomes, 264 also report OPEN unexcluded fruit/seed outcome and 167 have COMPLETE exclusion. Exactly **153 source experimental outcome quartets** contain DAY, NIGHT, OPEN and COMPLETE exclusion for the identical true published plant/site/year/season/metric units, from **79 independent studies / 85 plant taxa** (98 fruit-set + 55 seed-set). Source original HAND pollination appears in 126 matched day/night groups, of which 111 also have OPEN treatment; hand treatments are not automatically comparable standardized pollen supplementation assays.

## Four-arm functional definitions

Let O = open to both daytime/nighttime visitor access, D = day-only, N = night-only, C = complete visitor exclusion, all values original observed fruit or seed outcome.

**I = (O + C − D − N) / (O + C + D + N)**, bounded to [−1,+1] for finite nonnegative original source outcomes. Negative I = **subadditivity of original reproductive output** on the original measurement scale, not a direct behavioral demonstration of functional redundancy. Plant fruit/seed capacity saturation, abiotic pollen vectors and bag microclimate may all contribute.

**J = (O − max(D,N))/(O + max(D,N))**, bounded [−1,+1] when defined. Positive J = unrestricted day+night source outcome numerically exceeds the best measured one-window treatment. This is not proof that a specifically documented day visitor compensates for a moth loss.

Prior literature positive control **K = (O − C)/(O+C)**, which the original paper already interprets as pollination dependency; not a novel intervention finding.

## Actual original 79-study aggregate observations

To avoid 153 source site/year rows masquerading as 153 independent studies, the primary summary is the **median within each ORIGINAL published study**, followed by a sign test across the 79 study medians and an original-study-cluster bootstrap 3,000 draws. All intervals below are **descriptive** in the literature-selected studies and all p-values nominal unless corrected.

| Source original field reproductive endpoint | I nonadditivity | J open vs best | K open vs closed |
|---|---:|---:|---:|
| Original authentic four-arm comparisons | 153 | 153 | 153 |
| Independent original studies | 79 | 79 | 79 |
| Independent original plants | 85 | 85 | 85 |
| Median original study effect | **−0.05747** | +0.01450 | +0.83794 |
| Negative study medians | **49** | 31 | 4 |
| Positive study medians | 30 | 46 | 75 |
| Zero study medians | 0 | 2 | 0 |
| Exact study-level two-sided sign p | **0.04217** | 0.11002 | 5.24e−18 |
| Study bootstrap descriptive 95% interval of median | **[−0.09132,−0.01227]** | [−0.00106,+0.04706] | [+0.67624,+0.96010] |

**Multiplicity:** The two NEW primary exploratory endpoints are I and J. Holm adjustment of their genuine original-study p-values yields **adjusted I p = 0.08433** and **adjusted J p = 0.11002**. K is a known original-literature positive control, not a new discovery. Thus the recurrent negative I study direction is **not confirmed at multiplicity-adjusted 0.05**.

## Family-level and leave-whole-source sensitivity

Of 153 authentic four-arm outcome pairs, **127** have independently source-confirmed original botanical family membership from the original author 142-species botanical trait table, covering **29 families**. Only **9** families contribute three or more independently published studies; the median family has **two** study papers. The largest family has **26/127** original source field comparison rows.

- Across source family medians, **22/29 negative I**, **7 positive**, family median **−0.03143**, family-level nominal sign p **0.00813** (exploratory, family histories and sampling not independent randomized replicates).
- Remove each entire family: 29/29 overlapping original-study median I point estimates stay negative, range **−0.06336 to −0.04205**; not 29 independent experimental replications.
- Leave out each of **79 original studies** in turn: 79/79 source study-median I signs negative, coefficient median range **−0.06042 to −0.05465**, but only **30/79** overlapping source omission sign tests have nominal p<.05.
- Leave out each of **85 source plant species**: 85/85 original study-median I still negative, range **−0.06336 to −0.04848**, with nominal p<.05 in 44 of those 85 highly overlapping refits.
- Original botanical blooming-time support in known-family quartets: only **35 day** and **63 night** bloom study×site outcomes, with just **four original published papers** comparing both blooming categories. No original within-paper causal bloom mechanism established.

## Protocol, life history and actual study quality

| Independently specified real study restriction | Pairs / source studies | Study median I | Study sign p | Study bootstrap I median range |
|---|---|---:|---:|---|
| All 4-arm source experiments | 153 / 79 | −0.05747 | 0.0422 | [−0.0913,−0.0123] |
| Original reported native plants | 133 / 69 | **−0.03226** | **0.1480** | **[−0.1069,+0.0063]** |
| Original reported day/night exclusion time available | 88 / 52 | −0.06968 | 0.0365 | [−0.1069,−0.0177] |
| Fruit set outcome | 98 / 64 | −0.04669 | 0.0328 | [−0.1111,−0.0125] |
| Seed set outcome | 55 / 41 | −0.04698 | 0.0275 | [−0.0778,−0.0119] |

The native-only uncertainty interval **crosses zero**. The other post-hoc nominal subgroup p-values are correlated via common published studies; they cannot collectively prove a universal mechanism.

## Observational artificial nightlight: no supported difference

Original study-equal weighted WLS adjusted for experiment start year, geographic polynomial, original elevation and metric (fruit or seed), clustered by the 79 original published studies:

| Outcome × previous true Li–Cao satellite nightlight | beta fraction per prior mean log(1+PCNL) | Study-cluster CI95 | nominal p |
|---|---:|---:|---:|
| Four-arm I subadditivity | +0.0011203 | [−0.017495,+0.019736] | **0.9061** |
| J open versus better single period | −0.0025496 | [−0.025186,+0.020087] | **0.8253** |

Nothing here establishes that illumination alters plant reproduction, moth pollination, effective daytime replacement, selfing, or island adaptation. Original satellite nighttime light contrasts are largely between different papers and locations, not experimental lighting.

## Scientific conclusion and next source gate

The original experiment database provides a previously unquantified **cross-study four-arm source subadditivity pattern** in actual fruit and seed outcomes: 49/79 original journal studies, 22/29 families with source traits, but source-native restriction and Holm correction render a universal negative biological mechanism unconfirmed. Reproductive outcome can be subadditive through physiological fruit-capacity saturation, nonmatching bag effects or autonomous/wind background, not exclusively ecological functional redundancy. Also, unrestricted access does NOT robustly exceed the better single-window pollination treatment (J sign p=0.110).

The logical next experiment-level mechanism is to acquire actual **pollinator identity and effective pollen delivery at the same field site/year** from the independent 92 original papers flagged with author-reported visitor observations (PR #57, 91 DOIs). Species-level GloBI v0.11 positive guild records were found in just 23 of 123 source plants and NO original qualified positive relation has both event time and coordinates; those global records cannot resolve the missing same-site daytime/nocturnal guild comparison.

## Reproducibility and data stewardship

New source code: scripts/analysis/alan_original_diel_fourarm_reproductive_function_20261114.py; synthetic guardrails: tests/test_alan_original_diel_fourarm_reproductive_function_20261114.py; original SHA rehydration workflow: .github/workflows/alan-original301-diel-fourarm-reproductive-insurance-20261114.yml. Source study original receipt #43379308366, all 30 authentic PCNL raster annual source shards #43380089481, original author study Git blob verified, original main full-source beta unchanged. **Final SUCCESS run #43456561380**, artifact **#14054815172**, Draft PR #58. No original botanical study×location/year or source experimental raw treatment records exported; parent PRs and main unchanged.
