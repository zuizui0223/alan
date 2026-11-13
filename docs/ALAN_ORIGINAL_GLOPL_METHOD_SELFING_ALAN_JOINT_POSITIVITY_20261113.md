# ALAN: GloPL Sup versus Bagout × autonomous selfing × historical-light common-support audit

**2026-11-13. Original original-data replay SUCCESS.**  
[Full authenticated Github Actions run #43373289271](https://github.com/zuizui0223/alan/actions/runs/43373289271) · [aggregate-only artifact #14013930772](https://github.com/zuizui0223/alan/actions/runs/43373289271/artifacts/14013930772) · [Draft PR #50](https://github.com/zuizui0223/alan/pull/50).

## Purpose

Previous source-frozen GloPL 725 experiment analysis obtained a positive contrast **Sup − Bagout** for the interaction of preceding historical nightlight with independent species-direct autonomous selfing, **β=+0.86994 [+.37463,+1.36525], nominal study-cluster p=.00058**. Here we scrutinize **whether it is jointly supported within independent original study, same plant and true experimental year, and shared historical light**, rather than from nonrandom research design heterogeneity.

The original published [Bennett et al. 2018 GloPL data descriptor](https://doi.org/10.1038/sdata.2018.249) defines **Sup** as pollen addition to unbagged flowers, allowing natural pollination; **Bagout** as hand-outcrossing of bagged flowers, where natural pollination is excluded. Bagging may itself influence plant reproduction (e.g. microclimate). Original methods and geography are not randomized across the database.

## Original source identity and biological qualifications

- Source audit uses the original 1,419 actual pollen limitation experiments with real 1993–2015 experimental years (not article publication years), original 1992–experiment-year-minus-one MD5-authenticated Li–Cao annual PCNL, original photo-free GloPL reproductive log-response values. Source SHA of fixed prior-light join verified.
- Independent 2026-09-08 Island species-DIRECT high/medium and high-only autonomous versus absent classifications, original ledger SHA **7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577**. Unknown and multistate not assigned false biological categories.
- Source baseline source-wide actual-pollen-limitation versus historical light coefficient **β+.01071579263** reverified as checksum of same original published result. **The high+medium overlap is 725 original experiments, 242 accepted plant species, 179 original publications**, high-only 688/227/170. This is a mainland AND island GloPL sample, not an island-only source or actual moth-pollination trial.

## 1. Near-complete method/research confounding

| Original record relationship, high+medium source | Independent count |
|---|---:|
| Original published studies represented | **179** |
| Studies with only one method (Sup or Bagout) | **176** |
| Studies with BOTH methods | **3** |
| Studies with BOTH autonomous and non-autonomous plant species | **10** |
| Studies with BOTH methods AND BOTH autonomy classes (not necessarily each paired) | **2** |
| Botanical species that appear under BOTH methods anywhere | 3 |
| SAME original study×plant combinations using BOTH methods | **1** |
| SAME study×plant×TRUE experimental-year under BOTH methods | **0** |
| SAME study×plant×year×original site under BOTH methods | **0** |
| Studies with both autonomy classes and source PCNL range ≥0.03 | **4**, containing **131 original experiments** |

The lone source study×plant two-method pairing disappears under high-only botanical trait evidence. Direct within-study experiment-year comparison is **not possible**. Original treatment-method × trait × historical-light regression can still be fitted as a between-study **descriptive** relationship, but cannot establish a causal effect of bagging or pollinator exclusion.

## 2. Geographic and historical-light positivity failure

A 4-way source comparison requires measured original field experiment data for all combinations of **Sup versus Bagout** and **selfing yes versus no** within one genuinely comparable original geography AND some overlapping historical-light range. Merely including a geographic fixed effect does not create absent original observations.

| 5°, 10°, 20°, 30° region grid | 5° | 10° | 20° | 30° |
|---|---:|---:|---:|---:|
| Original regions | 110 | **75** | 39 | 27 |
| Regions containing all 4 source treatment×autonomy groups | 1 | **1** | 3 | 4 |
| Source cases in those four-class regions | 63 | **79** | 169 | 229 |
| Regions also having shared prior-PCNL exposure interval among all four groups | **0** | **0** | **1** | **1** |

At 10° spatial resolution the sole complete four-way block contains 79 rows, but the 4 groups have **no overlapping historic nightlight interval**. Even at a coarse 20° or 30° resolution, only one joint treatment×autonomy stratum has shared exposure support, far from geographically replicated factorial identification.

**Adjudication:** Old β=.86994 (p=.00058) is verified as a REAL ASSOCIATION of source methodology/plant/geography, **NOT evidence that the Bagout treatment changes how light and selfing alter field pollen receipt**. This distinction is not a p-value correction or a zero effect test; the needed common-support data do not exist in this source.

## 3. Actual nonbagged Sup pollen supplementation (measured reproductive output)

Primary source response: GloPL log reproductive response ratio of assisted pollen supply versus natural pollination, **not measured nocturnal moth-only pollination**. Positive values imply hand-added pollen increased reproduction. Model: study-equal WLS, original experiment year and geographic polynomial, prior-light exposure, measurement controls, **10° geographic FE**, published study-cluster nominal 95%CI; plant selfing categorical source is species-external, not measured at original flower location.

| Botanical-direct evidence and source outlier restriction | Experiments / studies / plants | Autonomous-minus-absent GloPL log-response difference β | 95% CI | nominal p |
|---|---|---:|---|---:|
| High+medium; Sup | **579 / 116 / 147** | **−0.18474** | [−0.44036,+0.07089] | .15664 |
| High+medium; Sup with abs(PL log ratio) ≤2 | **562 / 114 / 144** | **−0.19925** | [−0.40023,+0.00173] | .05200 |
| HIGH-only species-direct; Sup | **553 / 111 / 139** | **−0.09662** | [−0.35718,+0.16395] | .46737 |
| HIGH-only; Sup with abs(PL log ratio) ≤2 | **536 / 109 / 136** | **−0.09448** | [−0.28917,+0.10021] | .34151 |

Even the main Sup-only direct autonomous-selfing **main effect** is not statistically resolved under the stronger original high-quality evidence restriction; it cannot establish actual causal reproductive assurance. The earlier Sup-only **historical PCNL × selfing interaction** in PR #49 remained uncertain (β+.0453 [−.0705,+.1611] p=.443 high+medium; high-only β+.0038 p=.951). There is no consistently replicated nightlight-dependent selfing buffer among genuine measured Sup reproductive outcomes.

## 4. Source-method imbalance and sample support

| Original method×selfing source evidence | Bagout/no selfing | Bagout/autonomous | Sup/no selfing | Sup/autonomous |
|---|---:|---:|---:|---:|
| Actual experiment rows | **108** | **38** | **427** | **152** |

The original treatment methods occur in **66 Bagout source studies** and **116 Sup source studies**, with only three overlapping independent study identifiers. Bagout has **16/146** original study-case PL effects |log ratio|>2 versus Sup **17/579**. Original prior-PCNL mean was zero for 81/146 Bagout cases and 345/579 Sup cases. Source geographic/method support not random.

## 5. Scientific interpretation and next work

**Scientifically established:** genuine experimental pollen limitation differs descriptively by methodological protocol in the combined original data; joint support for method×selfing×historical-light mechanism is extremely sparse. Sup-only high-quality botanical records do not provide a stable significant effect of autonomous selfing or a stable historical-light interaction. The previously significant method contrast must be kept as a **source-design association, not a causal GloPL or moth-pollination effect**.

**Not established:** ALAN increases wild plant pollen limitation; Bagout itself causes source interaction; autonomous selfing compensates for loss of effective nighttime moth pollination; island night-pollination functional collapse or real plant-population shifts.

To make a direct mechanistic global analysis, the next source tier needs *same study and location*: before/after or paired illumination, pollen supplementation with consistent no-bagging rules, true day/night pollen transfer and mating system measured in the same wild population, including replicated sites and independent studies. More photographed site data or conditional geographic fixed effects alone cannot create missing experiment design contrasts.

## Reproducibility and privacy

Source-reread script: scripts/analysis/alan_GloPL_original_method_autonomy_common_support_20261113.py. Synthetic guardrails in tests/test_alan_GloPL_original_method_autonomy_common_support_20261113.py. GitHub Actions source SHA/dates, full replay and aggregate-only output: .github/workflows/alan-glopl-method-common-support-20261113.yml, run #43373289271. Draft PR #50, no mutation to parent PR #49 or prior ALAN science, and no publication or individual study×plant×field coordinates or rights-limited species trait values distributed.
