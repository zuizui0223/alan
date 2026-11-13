# ALAN: Direct experimental pollen limitation × previous nightlight × autonomous selfing

**Genuine original source replay completed on 2026-11-13, exploratory and noncausal.** [Final GitHub Actions run #43371281479](https://github.com/zuizui0223/alan/actions/runs/43371281479) SUCCESS; [aggregate JSON artifact #14011809833](https://github.com/zuizui0223/alan/actions/runs/43371281479/artifacts/14011809833). This is a new separate draft PR #49; existing PRs and main remain intact.

## Original independent reproductive experiment and historical light sources

- Source verified 1,419 original GloPL pollen-supplementation experiments at 597 field sites, 507 accepted plant taxa and 341 original studies, with real experimental years 1993–2015 (not publication year), and original Li–Cao PCNL nightlight from 1992 up through the year BEFORE each experiment. The original source-wide model's light slope was exactly reproduced: beta +0.0107158, CI [−0.0468274,+0.0682590], p=.715.
- Independently original signed 2026-09 Island high/medium directly documented species-level autonomous selfing, 134,268 botanical source trait ledger rows, 62,423 species across all trait dimensions, SHA256 7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577.
- Exact accepted plant-name and explicit autonomous / absent binary selfing only; unknown, multistate, facilitated and inferred genus labels are NOT coded as a negative.
- Actual independently qualified **725 pollen-limitation experiments / 242 accepted botanical species / 179 original studies**, with **65 autonomous species / 190 experimental cases** and **177 non-autonomous species / 535 experimental cases**. High-only botanically confirmed variant: 688 original cases, 227 species, 170 studies.
- GloPL responses are genuinely measured natural-versus-hand-pollination log reproductive output ratios, NOT local effective moth pollination. Source locations are global mainland and islands, **NOT an island-only population**.

## PCNL × direct autonomous selfing effect on measured pollen limitation

All models control for the original study source year, latitude and longitude polynomial, four GloPL measurement fields, original study-equal weights and study-cluster intervals. 10-degree geographic fixed effects were already a source baseline sensitivity, but the selfing moderation was studied posthoc.

| Study-provenance condition | n original GloPL experiments | Moderation beta | nominal 95%CI | p |
|---|---:|---:|---|---:|
| **Original broad geography/method/year controls** | **725** | **−0.13353** | [−0.29269,+0.02563] | **.1001** |
| 10-degree region fixed effects | 725 | −0.27901 | [−0.48955,−0.06846] | .00940 |
| 20-degree region fixed effects | 725 | −0.20673 | [−0.40322,−0.01025] | .0392 |
| **High-only direct botanical selfing, 10-degree** | **688** | −0.36662 | [−0.60909,−0.12416] | .00304 |
| Exclude extreme absolute pollen-limitation logratio >2, high+medium | 692 | −0.08950 | [−0.20065,+0.02165] | .1145 |
| High-only exclude same extremes | 655 | −0.13160 | [−0.24143,−0.02177] | .01886 |
| **Normal nonbagged pollen supplementation (Sup) only** | **579** | **+0.04531** | [−0.07045,+0.16107] | **.4430** |
| High-only Sup only | 553 | +0.00380 | [−0.11753,+0.12514] | .9510 |
| Original sites with nonzero prior PCNL only | 299 | −0.06859 | [−0.66940,+0.53222] | .823 |
| 10-degree blocks plus botanical genus FE | 725 | **NOT ESTIMABLE (rank deficiency)** | — | — |

Study-and-botanical-species TWO-WAY cluster sensitivity for the 10-degree 725-case moderation remains beta −0.2790, finite-t CI [−0.49092,−0.06709], p=.01016. The high-only two-way test yields p=.00350. The original primary is still NOT statistically supported, and Sup-only remains non-significant. Do not promote the conditional 10-degree finding as evidence of ecological ALAN damage/rescue.

**Source-independent overlap challenge:** only 10 of 179 original GloPL papers have both selfing categories. 75 original 10-degree grid blocks include plants; 24 have both categories; only 19 also have genuine within-block light standard deviation >.03 (351 experiments, 87 studies). Some 58.8% of the 725 source experiments have a historical PCNL mean exactly zero.

## Full one-at-a-time published-study and plant influence

Under the original 725-case 10-degree model only: all **179/179 published study omissions** and all **242/242 plant species omissions** yield a negative point estimate; both ranges [−0.31611,−0.21932], median approximately −0.2790. Thus no one published study/plant causes the conditional point sign. These highly overlapping deletions are **not independent replication** and cannot repair source method or genotype/region confounding.

## Critical NEW source-method heterogeneity

The original GloPL authors explicitly note TWO biologically different hand-pollination protocols:
- Sup = additional pollen given to unbagged flowers, natural pollen arrival still possible;
- Bagout = bagged, artificially outcrossed flowers, excluding natural pollinators during treatment. Bagging can change seed output independently through flower microclimate.

See the original GloPL data descriptor, https://www.nature.com/articles/sdata2018249 . The GloPL source field PL_Effect_Size_Type2 distinguishes these original protocols.

| Source hand-pollination procedure | Genuine cases | Accepted species | Studies | Autonomous/not-autonomous source cases |
|---|---:|---:|---:|---:|
| Bagout | **146** | 98 | 66 | 38 / 108 |
| Sup | **579** | 147 | 116 | 152 / 427 |

A simultaneous **ALAN × documented autonomous selfing × actual GloPL METHOD TYPE** diagnostic is fitted with common source geography/year/method controls. The reported coefficient compares the PCNL-autonomy moderation SLOPE under Sup MINUS that under Bagout, not the absolute selfing effect or the GloPL source baseline:

| Source-provenance and alternative original method contrast | beta Sup minus Bagout | nominal 95%CI | p |
|---|---:|---|---:|
| **High+medium direct selfing, 10-degree blocks** | **+0.86994** | **[+0.37463,+1.36525]** | **.000577** |
| High+medium 20-degree blocks | +0.72758 | [+0.28070,+1.17446] | .00142 |
| High+medium excluding PL logratio abs>2 | +0.40952 | [+0.15304,+0.66600] | .00175 |
| High-only direct trait, 10-degree blocks | +0.93093 | [+0.39686,+1.46500] | .000635 |
| High-only excluding PL logratio abs>2 | +0.46329 | [+0.21828,+0.70830] | .000211 |

Original study and species two-way clustered p for the first difference=.000723. **This exploratory source-design heterogeneity is more consistent than the claimed ALAN ecological effect**. It does NOT prove bagging itself caused the slope: the Bagout and Sup experimental sources are not randomized within identical original study/population samples, and experimental methods are correlated with publications and geography. This is a methodological identification problem, not experimentally demonstrated physiological response to light.

## Independent nocturnal functional evidence is too sparse

Original independent published day/night reproductive-exclusion studies identify 134 plant species; only **11 accepted plants** match GloPL's measured pollen limitation AND the strict botanical direct autonomous-selfing flag: **64 GloPL cases / 16 original studies**. Of these 11, all have directly documented autonomous selfing ABSENT, while 7 have numerically day-higher versus 4 night-higher reproductive responses in their source studies. Therefore a true **prior PCNL × actual night function × autogamy** comparison is **NOT ESTIMABLE**. The independent night-pollination exclusion experiment cannot be assumed to identify moth-only pollination.

## Adjudication and data needs

This analysis establishes (a) the strong dependence of historical-light/autogamy moderation estimates on spatial and GloPL treatment definitions, and (b) measurable source method heterogeneity of study comparisons. It does **not** establish more pollen limitation in night-light-exposed nonselfing plants, autonomous seed-set rescue, effective moth-specific pollination loss, nor an island-specific interaction. Pollen limitation is genuinely measured, but the direct autonomous-selfing trait was reported for a species outside the GloPL field population; PCNL is remotely sensed previous-year light, not local manipulated illumination; most independent studies only measure one selfing type.

Priority for biological mechanism: geographically and temporally paired original natural-versus-supplemented pollen experiments, effective day/night pollen delivery and reproductive outcomes, autonomous selfing measured at same wild plant population, and independently measured lighting or replicated restoration. Pre-register focal methods/replicates rather than p-value selection.

## Reproduction

- [Signed source actual run #43371281479](https://github.com/zuizui0223/alan/actions/runs/43371281479)
- [Final aggregate-only JSON artifact #14011809833](https://github.com/zuizui0223/alan/actions/runs/43371281479/artifacts/14011809833)
- Source pipeline scripts/analysis/alan_glopl_direct_autonomous_selfing_PCNLPrior_20261113.py; original 1419 source baseline helper scripts/analysis/alan_glopl_real_PL_priorPCNL_exploratory_20261110.py; tests/test_alan_glopl_direct_autonomous_selfing_PCNLPrior_20261113.py; workflow .github/workflows/alan-glopl-autonomous-selfing-priorPCNL-original-20261113.yml. Full SHA gates and no source experiment plant×location or rights-limited botanical trait rows re-published. Earlier PRs and main unchanged.
