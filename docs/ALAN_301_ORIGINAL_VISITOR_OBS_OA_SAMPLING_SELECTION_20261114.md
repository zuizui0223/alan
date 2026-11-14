# ALAN original 301 real reproductive outcomes — observation reporting and OA access selection audit

**Verified original-source run SUCCESS**: [GitHub Actions #43561645268](https://github.com/zuizui0223/alan/actions/runs/43561645268), [aggregate-only artifact #14100451346](https://github.com/zuizui0223/alan/actions/runs/43561645268/artifacts/14100451346), Draft PR #62. Both original Wenzell direct 2019 population-year source and Marques original bat/bird published fruit/seed source table were reverified in the same successful run.

## Experimental source and study selection

Original 301 true day/night fruit/seed reproductive comparisons arise from 118 original studies / 123 plant species; signed pre-PCNL original source SHA256 21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af and original 1,652 treatment published author Git blob 39a82d966afcdde670603add1c881ace8791ffe9.

The authors report some pollinator observations in **92** of those original studies (247 of the 301 field experimental outcome comparisons). The other **26** studies (54 outcomes) have no positive source observation flag or unknown flag; **not** proof of local biological visitor absence. Separately, among the 92 flag-positive studies, OpenAlex's independently retrieved bibliographic metadata identifies **36 candidate OA links**, corresponding to 105 original reproductive outcomes, versus **56 without available OA candidate**, containing 142 original reproductive outcomes. OA links indicate access to article text, **not** actual local moth or butterfly events.

## Original real observed night reproductive fraction by source selection

Response is the original observed night/(day+night) fruit-or-seed response, first summarized by median within each **independent published study**, then averaged across original source papers. Each study receives equal weight, not each field plot or photo; 3,000 publication-level bootstrap resamples produce descriptive intervals. Groups are pre-defined by original author observer flag and OA public metadata; no source outcomes or p values are used in selecting cohorts.

| Study cohort | Original independent papers | Measured original reproductive cases | Mean of study-median night fraction |
|---|---:|---:|---:|
| Positive source observer flag | 92 | 247 | 0.57501 |
| Other / unknown observer flag | 26 | 54 | 0.55895 |
| Observer flag YES + OA link candidate | 36 | 105 | 0.61386 |
| Observer flag YES + no OA candidate | 56 | 142 | 0.55003 |

**Observer flag 92 versus 26:** source night-share difference +0.01606, publication-block bootstrap descriptive 95% **[−0.10183,+0.13851]**. No robust night reproductive source cohort shift established. Original median study experiment years 2008.5 versus 2010.25; absolute latitude medians **25.88° versus 36.15°**.

**OA 36 versus other observer-flag 56:** source night-share difference +0.06382, publication-block descriptive 95% **[−0.02874,+0.16082]** (also includes zero). Original median experiment years **2011.5 versus 2005.5** and median absolute latitude **30.31° versus 24.95°**: observable sample composition differs, although no causal OA effect has been estimated.

**Consequence:** Taking only 36 easily accessed OA papers drops **142 of 247** actual field source reproductive comparisons in the original 92 observer-flag studies. That OA subset is not automatically globally representative for diel pollinator function or historical artificial-light exposure. Use source publication/year/geography missingness gates, prioritize collecting real guilds from non-OA source papers where feasible, and preserve within-study biological independence rather than maximize a positive coefficient.

## Direct verified ecology context

**Wenzell 2024:** one original independent publication, with real source 2019 visitor-dominant guild and plant-level fruit set in **7** observed population-years (204 original plant treatment records). **6** of these are in the authentic 301 meta-analysis, 5 strict spatial matches within 1 km, 1 provisional spatial match 1.65 km. Published original ecological contrast in all 7: mean nocturnal night/(night+day) fruit share **0.69845** (3 hawkmoth-dominated populations) versus diurnal **0.24941** (4 bee/hummingbird-dominated populations), +0.44904 descriptive, **not seven independent studies**.

**Marques 2015:** a second independent original paper, with bat/hummingbird visitor identity recorded in the same field research area and reproduced original primary published fruit/seed experimental source table, night-only fruit set **647/788 = 82.1%** versus day-only **81/808 = 10%**. Exact visitor-event year not independently confirmed against original 2010 global source experiment label; this is *not* another fully year-matched visitor study.

The original global 301 study historical prior-PCNL night-share model remains unresolved (β +0.00307, 95%CI −.02239 to +.02852; p=.813). No randomized ALAN experiment, ecological moth extinction, daytime compensation, or island reproductive response established.

## Reproduce and stewardship

Source-neutral reproducible code: scripts/analysis/alan_original301_observer_OA_sampling_bias_20261114.py and associated synthetic unit tests. Original independent public 92-paper bibliography run #43538019810 authenticated. Full original source job #43561645268 SUCCESS with aggregate data artifact #14100451346. No original study×plant×site/year rows or original animal observation-event coordinates exported. Historical main and old PRs untouched.
