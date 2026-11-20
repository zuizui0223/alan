# Grognuz et al. (2026) 14-paired-site ALAN seasonal phenology: can a wet-plant trait moderator be tested?

2026-11-20 | ORIGINAL PAPER + INDEPENDENT UK FLORA source verified | DATA SOURCE/TRAIT COMMON-SUPPORT STOP | \`main\` unchanged

## Primary objective

Keep the existing long-term 3,901-plant candidate distribution result ("species classified as lit-associated recording decliners have higher wet-habitat moisture affinity"; adjusted +0.534 EIV points) as the paper's main **observational** finding. Investigate why ALAN might affect wet-associated taxa by independently testing controlled lighting experiments, but do not claim a causal pathway simply because lighting changes flowering in some other plants.

This audit asks whether [Grognuz, Gisler & Knop (2026), *Biological Conservation* 314:111650](https://doi.org/10.1016/j.biocon.2025.111650), an externally controlled night-light experiment, provides enough wet vs drier plants to estimate a **moisture × experimental ALAN flowering timing** interaction.

## Grounded original study design and endpoint

[Official Agroscope archived author article PDF](https://ira.agroscope.ch/en-US/Page/Einzelpublikation/Download?einzelpublikationId=68011), original Fig. 4 labels and methods:

- 2022 = 6 lit/dark site pairs; 2023 = 8 pairs, nominal total **14 independent paired sites over two years**, NOT 108 independently assigned light treatments.
- For **seasonal flower phenology**, in 2022 one dark site was excluded: 11 instead of 12 sites. 2023 used 16/16. Real flowering experiment contained **108 source phenology plots** (4 plots per observed site), with plant taxon × plot repeated weekly. Source site-pair match keys are required before treating all 14 nominal pairs as 14 complete flower-phenology contrasts.
- Sown mixtures of 22 annual wildflowers on plots. **16 plant labels** had ≥10 lit and ≥10 dark plots observed and were kept by authors for phenotype modelling, original paper Fig. 4.
- **Four source endpoints are timing** in days after beginning of seasonal spring: FIRST, QUARTER (one-quarter peak), PEAK, and ULTIMATE (last open inflorescence). They are **not** maximum flower count, daily flower surface/cover, compatible pollen delivery or seed success. Authors found ≥1 timing shift in 12/16 species, with both advances and delays; this is their **published** result, not a new dataset computed here.
- Original author Mendeley deposited dataset is labelled [10.17632/xjmn5jhxyn.2](https://doi.org/10.17632/xjmn5jhxyn.2) in the paper (DOI version 2 described as reserved); the [public Mendeley dataset page](https://data.mendeley.com/datasets/xjmn5jhxyn/1) currently exposes version 1, 10.17632/xjmn5jhxyn.1.

## New original Fig. 4 × independently maintained Kew plant-moisture trait audit

[Successful first original data+source gate CI #44673329871](https://github.com/zuizui0223/alan/actions/runs/44673329871), 3 synthetic source tests PASS, original figure label list validated against author article. [Saved exact 16-source-name, F values and API statuses](../results/all_species_forest/biotic_interactions/grognuz2026_16plant_Fig4_KewF_Mendeley_gate_20261120.json).

RBG Kew [BIFloraExplorer *raw* British–Irish vascular plant trait CSV](https://github.com/RBGKew/BIFloraExplorer/blob/master/data-raw/BI_MAIN.csv), authentic source Git blob SHA1 \`e6a1da929af81f370d2a1b79177b489e3e573a06\`, **3,227 plant-taxon rows**, exact original \`taxon_name_binom\` only. UK moisture \`F_PLANTATT\` versus independently compiled \`F_Doring\` retained **separately**, not pooled or silently filled. An F score is a regional ecological preference indicator, NOT a measured local Swiss water table.

| Original study Fig. 4/UK moisture support | Verified source count |
|---|---:|
| Original flower-phenology taxon labels | **16** |
| Exact original Kew binomial matched | **13** |
| UK PLANTATT F present | **12** |
| UK Döring F present | **10** |
| UK PLANTATT F >7 (wet species among source-scored) | **0** |
| UK Döring F >7 (wet species among source-scored) | **0** |
| UK PLANTATT F range | **2–5** |
| UK Döring F observed values | **2,4,5** |

All 12 truly scored original UK plant taxa are **nonwet F≤5 on the PLANTATT indicator**. Four unscored original labels are **NOT classified as nonwet**:
- \`Sonchus sp.\` and \`Epilobium sp.\`: unidentified below genus; no exact species ecology can be assigned.
- \`Silene pratensis\`: no exact binomial record in the pinned Kew source; do NOT substitute \`Silene latifolia\` without authoritative synonym data.
- \`Centaurea jacea\`: exact Kew binomial exists but UK F score is \`NA\`.

This is a strong **design-based negative source-support result**: the original study can establish ALAN-driven phenology changes for the source *wildflower-strip assemblage*, and may eventually test **variation within UK F=2–5**, but it simply does **not** contain a UK-defined wet group F>7 with which to estimate wet-versus-drier experimental sensitivity. It cannot validate or refute the main wet-species macro bias.

Neither 108 plot measurements nor 16 independently modelled plant response coefficients can create four wet species that are absent from the scored original source. Formal \`light × F\` or \`light × wet\` phenology coefficients are **NOT fitted**.

### 同じ原著の発芽・花弁の日周運動まで形質を照合

[追加原本検証CI #44673329871 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/44673329871)で、著者論文中の別々の実験ステージも英国Kewの同じ真正原CSVから独立に検索した。

| 元研究の別々の応答 | 原論文の植物種 | PLANTATT水分値を確実に照合できた種 | F>7の湿潤種 |
|---|---:|---:|---:|
| 開花フェノロジー（FIRST/QUARTER/PEAK/ULTIMATE） | **16** | **12** | **0** |
| 発芽タイミング（別のトレー／種子実験） | **6** | **4** | **0** |
| 花弁の開閉時刻（日周運動・短期/長期光曝露） | **2** | **2** | **0** |

発芽で形質値があったのは *Cichorium intybus*（F4）、*Crepis capillaris*（F4）、*Hypochaeris radicata*（F4）、*Trifolium incarnatum*（F2）の4種。*Anthemis tinctoria* と *Silene pratensis* は真正Kew二名法値がないため乾性に割り当てない。花弁開閉の *Silene noctiflora*（F4）と *Malva moschata*（F3）は双方スコアを得たが、いずれも元の英国指標では湿潤性7超ではない。

**16+6+2=24独立植物種とはならない**。実験プロトコルは別々で、植物種の重複（例：*Malva moschata*）もあり、ランプ実験の物理単位は根本的に共有しうる。また原著の12/16開花タイミング変化を「湿潤植物が繁殖を落とした」ことと混同しない。

### Source Mendeley public API was actually tested

The [currently documented official public API](https://data.mendeley.com/api/docs/) at \`api.data.mendeley.com\` was queried via GitHub Actions on 2026-11-20:
- \`/datasets/xjmn5jhxyn\` (latest), \`?version=1\` and \`?version=2\`;
- \`/datasets/publics/xjmn5jhxyn/files?version=1/2\`;
- \`/datasets/xjmn5jhxyn/versions\`.
- **All six official JSON data/file/versions endpoints returned HTTP 401** on this runner. The public data.mendeley.com HTML v1 page returned HTTP200 but **not raw data or a validated file manifest**.

Therefore zero source experimental plot rows, original site treatment-ID assignments, or within-species variance have been downloaded or fitted. We neither substitute a screenshot-digitized Fig. 4 point estimate nor pretend the figure's four phenology stages form independent light experiments. ACCESS STOP is separate from the plant-trait COMMON-SUPPORT STOP.

## Relation to prior independent mechanism evidence

| Original programme | Source replicate design | Independent Kew/European F score high-wet support | What is identified |
|---|---|---|---|
| Giavi 2021 original 6 lit/dark meadow pairs; 1,884 underlying observed insect group records | Only **2 original physical meadows** with simultaneous wet+drier named taxa in the source; 3 or 7 paired species-rich survey events depending on F dataset | Two scored wet plants | ALAN-time visitation/flower-cover associations are exploratory; original floral-cover association survives one-plant-out, but **dry-only F≤7 slope also negative**, no wet-specific effect |
| Anic et al 2022 13 meadow species, 54 real physical troughs, one lighting programme shared with older Bennie papers | Source 665 plant×trough records and original 3 light codes (NOT mapped to colors) | **Only one** F>7, *Lotus pedunculatus* (UK F=8) | Original species-light biomass differences; Kew F × light-code rank slope -0.634, Holm/max3 p=.0766, not confirmed |
| **Grognuz et al. 2026** 16 phenology species, total 14 physical 2022/23 lit/dark pairs with one 2022 dark phenology site missing | 108 nested seasonally observed source plots | **Zero** of 12 UK F-scored plants F>7, four source phenology labels unscored | Author-reported 12/16 shifts in *flowering timing*, not wet-specific ALAN |

This does not demonstrate that *no Swiss plants are wet-habitat plants*. The separate UK flora trait system may miss local habitat ecology, and four published taxa are not scored. It demonstrates that **these particular published experimental source labels and independent available UK scores do not provide identifiable wet-vs-drier common support**. No multi-experiment wetness×ALAN seed/flower or pollen outcome has been established.

## Manuscript decision

1. Retain primary discovery about over-representation of high-moisture-affinity macro recording decliners. It remains observational, selection-sensitive, and not proof of ALAN-caused plant decline.
2. Use Grognuz in the **experimentally observed altered timing of floral resource availability** mechanism branch, distinct from Giavi's contemporaneous *flower-cover area* results.
3. Do NOT put the Grognuz 16 species × 4 stage × 14 nominal site counts into a false 896-independent light-effect meta-analysis; do not fit nonexistent wet trait moderators.
4. More promising candidate experimental work requires site-level water or matched wet/drier plant responses in original illumination treatments (and direct flowering amount, compatible pollen or seed output). A within-species ALAN × **manipulated soil moisture** experiment (Hey et al. 2020 DOI 10.1007/s00442-020-04679-9) shows direct water interaction with **larger positive ALAN growth effects under no watering**, i.e., it does not straightforwardly predict wet-species declines, and is only one species. Live-macrophyte Xue et al. 2026 ALAN×grazing is a separate four-aquatic-species design with no dry comparison and independently inaccessible original biomass rows as yet.

**Hold as a Draft PR, main unchanged.** This is an inspectable exclusion criterion supported by authenticated species labels and trait source, not an invented negative light biology result.

## Reproduction
- [16 original Fig. 4 names, Kew trait and official Mendeley source probes](../scripts/acquisition/alan_grognuz2026_fig4_kewF_mendeley_source_gate.py)
- [Synthetic Kew 12-species / missing-F / HTTP401 guards](../tests/test_alan_grognuz2026_fig4_kewF_mendeley_gate.py)
- [Original source re-execution workflow](../.github/workflows/alan-grognuz2026-original16-moisture-20261120.yml)
- [Exact output JSON](../results/all_species_forest/biotic_interactions/grognuz2026_16plant_Fig4_KewF_Mendeley_gate_20261120.json)
