# ALAN — 夜間光の増加と湿潤な生育地の植物

## 主要結果（2026年11月17日固定）

増光した地点で記録率が低下した候補283種には、調査地点の気温・降水・市街化の変化をそろえても、湿潤な生育地を好む植物が多かった。1992–2024年の長期記録で評価した3,901種のうち、湿潤性指標が得られた候補は160種。調査地域、地点数、島の地点割合をそろえて比較できた158種では、平均湿潤性指標は候補6.09、参照群5.50（差＋0.59、17指標のBH補正後 q=0.030）だった。同じ158種について増光とこれらの環境変化を同時に調整した差は＋0.53（主要調査地域ジャックナイフによる近似95%区間＋0.16～＋0.91）。湿潤性指標が7を超える割合は35.4%対22.3%だった。長期調査地点は欧州に偏り、参照群を増光地点で減少した他種に限ると平均差は＋0.35（名目 p=0.083）だった。これを本研究で唯一の主要結果とする。[原稿](manuscript/ALAN_WET_HABITAT_ja.md)・[固定した283種](results/all_species_forest/decline_traits/cohort.csv)・[形質別の推定値](results/all_species_forest/decline_traits/traits.csv)

### 論文用の結果と考察

[論文草稿の要旨・結果・考察](manuscript/ALAN_WET_HABITAT_ja.md)に、湿潤性の群間差、ホソヒラタアブとカニグモを含む共通経路、蛾の体表花粉、内陸湿地の年次観測を集約した。湿潤な減少候補に含まれる *Cirsium oleraceum* では、独立した照明実験の再解析でも昼の訪花が減る方向を示したが、草地6ペア単位では不確かさが大きかった。[実験の再解析と反例](results/all_species_forest/wet_environment/README.md)

## 副次的な探索結果

**島と本土の比較（旧主要結果）**　2019・2021・2023・2024・2025年の78地域、10,698地点年で、観察された植物種に占める蛾の訪花記録を持つ種の割合について、夜間光に対する島−本土の傾き差は **−1.902** パーセントポイント／log(1＋VIIRS) だった（地域クラスタ名目95%信頼区間 −3.706～−0.098、未補正 p=0.039）。各種を地点年ごとに一度数えた群集構成の推定である。GloBI v0.11で種名と訪花記録を更新した同地点比較では負の差が維持されず、結果は使用する記録集合に依存した。[元の推定と固定条件](https://github.com/zuizui0223/alan/blob/ff069a160400588b8d6134c48b5481729d75eb9a/docs/ALAN_GLOBI_PRIMARY_RESULT_FREEZE_20261111.md)・[記録更新後の比較](https://github.com/zuizui0223/alan/blob/ff069a160400588b8d6134c48b5481729d75eb9a/docs/ALAN_GLOBI_V011_SAME_FIVEYEAR_20261111.md)

**花粉制限（副次結果）**　島と本土に共通する544種のうち、補助授粉実験がある64種では、花粉制限が強い種ほど、同じ78地域で夜間光に対する種別の島−本土の傾き差が負に大きくなる傾向があった。花粉制限指標（手授粉／自然授粉の繁殖成績の対数比の種別中央値）が1高いと、傾き差は **−0.0194** ポイント変化した（地域クラスタ95%信頼区間 −0.0382～−0.0005、未補正 p=0.044、2仮説のHolm補正後 p=0.088）。花粉制限は文献の実験値であり、撮影地点での送粉不足や照明の因果効果は未検証である。[推定値](https://github.com/zuizui0223/alan/blob/ff069a160400588b8d6134c48b5481729d75eb9a/results/mating_light_contrasts/contrasts.csv)・[生態学的考察](https://github.com/zuizui0223/alan/blob/ff069a160400588b8d6134c48b5481729d75eb9a/manuscript/ALAN_GLOBI_ISLAND_POLLINATOR_SHIFT_DISCUSSION_ja.md)

## 機構仮説

同じ158候補で環境変化を同時に調整したマクロ比較と、湿潤な減少候補のアザミを含む照明実験の昼間訪花データを再解析した。湿地水文、水面反射、食害・捕食、送粉を区別した一次研究の検討を含む。[追加の検証結果](results/all_species_forest/wet_environment/README.md)

水面反射、訪花・送粉の変化、捕食を介した波及は、湿潤種の偏りを説明しうる仮説として扱う。既存の照合では同じ地点年の捕食と訪花を結べず、これらの寄与は未推定である。[検証範囲と出典](manuscript/ALAN_WET_HABITAT_ja.md)

湿潤植物363種の送粉・繁殖様式を照合し、同一調査地点での昆虫媒・風媒の記録率傾きを比較した。蛾の花粉運搬を確認できた減少候補7種について捕食記録との接続も整理した。[種別の検証](results/all_species_forest/wet_reproduction/README.md)

## 旧研究の探索結果・解析履歴

以下に花色・島面積・写真資料に関する過去の解析と再現情報を保存する。

<details>
<summary>旧研究の問い、結果、再現手順を表示</summary>

### Historical floral-colour question

Does the spatial association of artificial light at night (ALAN; VIIRS upward radiance) with the *photographed flower-region white-colour fraction* differ between island and mainland sites?

### Historical floral-colour estimate

- Five observation years: **2019, 2021, 2023, 2024, 2025**.
- **7,982** source-backed 0.05° photo-cell × year records, including **710** island cell-years.
- Island × log(1+VIIRS) spatial interaction: **−2.497 pp** (nominal spatial-cluster 95% CI **−3.949 to −1.046**; p **0.000802**), after climate, elevation, population, photography effort, geographic and observation-year controls.
- On the broader observed-photo eligibility cohort (20,543 cell-years), the estimate attenuates to **−1.205 pp** (95% CI **−2.487 to +0.077**); matched same-year/geographic comparisons and repeated sites do not establish a universal island-specific light response.
- The measured endpoint is an **original photographic floral-ROI white component**, NOT independently verified white-flowered species prevalence, causal trait evolution, effective nocturnal pollination, or local extinction.

The original source runs, numerical results and caveats are preserved in [the migrated draft](manuscript/MOIS_ISLAND_ALAN_OBSERVATIONAL_COLOUR_ECOLOGY_20261108_RESULTS_DRAFT.md) and [analysis notes](docs/MOIS_ALAN_ISLAND_FIVEYEAR_EMPIRICAL_20261107.md).

## Structure and execution

- `manuscript/`: draft observational Results/Discussion.
- `docs/`: primary empirical and data-source audit notes.
- `scripts/`: original research scripts, **preserved under their existing import paths**.
- `tests/`: historical tests; successful original GitHub run links remain in the source notes.
- `results/`, `config/`: published aggregate tables, checksums and source queues; **not the original photo/VIIRS data archives**.
- `archive/mois_workflows/`: copied legacy GitHub Actions as **non-executing historical recipes**. They refer to MOIS action artifacts and old branch triggers; do not simply move them into `.github/workflows/` without adapting and validating external inputs.

## Migration / provenance

See [MIGRATION.md](MIGRATION.md). The original signed image/data artifacts remain in the source MOIS release and Actions history; source files have not been remeasured and missing colour values have not been filled. Scientific analysis continuation belongs here, and moisture–flower-colour analysis remains in MOIS.

## Native ALAN reproduction (verified)

On 2026-11-08, ALAN's own [five-year source-verified GitHub Actions run](https://github.com/zuizui0223/alan/actions/runs/42438862639) succeeded using the historical MOIS original-photo and same-year VIIRS/WorldPop source artifacts. It reproduced **7,982 photo cell-years (710 island)** and the island × log(1+VIIRS) white-ROI spatial association **−2.497045 pp**, 95% CI **[−3.948547,−1.045542]**. The stricter same-year same-geographic island/mainland comparison was **not individually resolved**; no light-caused biological white-flower loss is inferred.

Use `.github/workflows/alan-native-fiveyear-replay.yml` via the **Actions → Run workflow** manual control on the default branch. This workflow intentionally leaves the historical MOIS workflows inactive. See [full original-run outputs, safeguards, and archival deadlines](docs/ALAN_NATIVE_REPLAY_VERIFIED_20261108.md).

## Exploratory island-area result (2026-11-08)

A source-verified new ALAN run on real GSHHG area and original flower-region ROIs found a **positive area × nightlight interaction** for within-island photographic whiteness: **+2.209 pp per tenfold island area** (nominal island-cluster 95% CI **+0.007 to +4.411**, 44 spatially supported islands/650 photo cell-years) in the strict cohort, but **+0.998** (CI **−0.390 to +2.387**, 95 islands/1,884 cells) in the broader measured-photo cohort. The source-support audit found that the larger-area island third supplied ~84% of photographed cell-years in both cohorts. Giving **each island equal total weight** produced **+5.152 pp** (island-cluster nominal CI **+2.559 to +7.744**) in the strict 44-island cohort and **+3.385 pp** (**+1.319 to +5.450**) in the broad 95-island cohort. These estimates support a **positive observational area moderation** under both island-equal cohorts, but are exploratory and not independent causal proof. Thus a smaller-island stronger-negative-light slope is **exploratory and sensitive to cohort support**, not a universal ALAN effect or confirmed pollinator mechanism. See [original empirical results](docs/ALAN_GSHHG_ISLAND_AREA_LIGHT_EMPIRICAL_20261108.md) and [successful ALAN run](https://github.com/zuizui0223/alan/actions/runs/42439477013).

## Island area versus mainland distance (verified)

The new source-verified [coastline-isolation analysis](docs/ALAN_MAINLAND_ISOLATION_VS_ISLAND_AREA_EMPIRICAL_20261108.md) compares island size with nearest GSHHG mainland coastline distance, using **the same 2019–2025 original photographed floral ROI source**. In joint within-island, **island-equal** models, the island-area × nightlight interaction remains positive in both the strict **44-island** cohort (**+5.384 pp**, 95% CI **+2.678 to +8.089**) and broad **95-island** cohort (**+2.939 pp**, 95% CI **+0.833 to +5.045**). The mainland-distance interaction is uncertain in both. **Cell-year-equal area intervals also include zero**; this remains a sampling-sensitive, exploratory spatial association, not causal evidence for flower evolution or nocturnal pollination changes. [Reproducible ALAN run](https://github.com/zuizui0223/alan/actions/runs/42440394223) · [signed result summary](results/alan_isolation_area_light_20261108_source_verified_summary.json).

## Photographed taxon composition and island area × nightlight (new real result)

ALAN has now [source-verified and executed a same-island, same-photo-cell-year decomposition](docs/ALAN_ISLAND_AREA_TAXON_COMPONENTS_EMPIRICAL_20261108.md) for the positive island-area × nightlight observation after adjusting continental isolation. In the **equal-island** model, strict 44-island interaction **+5.384 pp = +4.296 pp** from observed accepted taxa's *same-photo-archive* white-reference component **+1.088 pp** from taxon-relative photographed colour; broad 95-island **+2.939 = +2.059 +0.880 pp**. These are exact algebraic photography-based splits, **not biological species replacement or causal mediation**.

A more independent **2019/2021/2023 photographic species-colour reference evaluated in 2024/2025** has support in **29 strict / 50 broad islands**, but the island-equal earlier-year taxon-reference component has intervals overlapping zero. Thus **species assembly cannot yet be identified as the unique ecological mechanism**. Within-taxon photographic deviations remain possible but include observer/season/ROI effects, not verified flower-colour evolution.

[ALAN execution #42442841980](https://github.com/zuizui0223/alan/actions/runs/42442841980) · [Compact source-signed result](results/alan_island_area_taxon_components_20261108_verified_summary.json). The earlier island-area and coast-isolation results remain unchanged.

## Real additional floral ROI information from already measured sibling photos (2026-11-08)

The signed MOIS 1,181,146-original-photo / 451,643-strict9 archive contains **148,264 extra photo IDs** beyond ALAN's canonical 1,032,882 original one-photo-per-observation reference. They are all **other photos of existing observation UUIDs, not new independent observations**. Of 49,940 genuine extra strict9 photographs, **18,950 existing observation UUIDs** have a strict9-successful sibling when their original representative photo was technically unknown.

An ALAN-native SHA/GSHHG/WCVP/source-year audit found that admitting these genuine sibling-photo colour results **without increasing photo-opportunity or observation counts** increases **strict-eligible photo cell-years 7,982 → 8,515 globally and 710 → 793 on islands (+83, spanning 25 islands)**; the 3-photo/2-accepted-species cohort rises **20,543 → 21,680 globally and 2,026 → 2,155 on islands (+129 across 41 islands)**. Forty geographical cells with inconsistent original accepted-taxon reconciliation were excluded from additional admission, not removed from the prior results. **No new island×nightlight regression was fitted on the expanded source:** the values are genuine measured-photo opportunity gains, not a new ecological effect.

[Complete real source evidence](docs/ALAN_MORE_AUTHENTIC_SIBLING_PHOTO_ROI_SOURCE_20261108.md) · [Successful ALAN data-only run](https://github.com/zuizui0223/alan/actions/runs/42445111291) · [aggregate JSON](results/alan_20261108_extra_measured_sibling_ROI9_photo_opportunities.json).

**To add truly independent photo observations next**, select new observation UUIDs from the already-indexed 93.2m-photo upstream candidate frame before downloading any image pixels; keep source-reference and targeted refill selection separate. The MOIS repository remains moisture-focused.

## Real acquisition from the full 93.2m photo index (2026-11-08)

The original-photo source is now expanded beyond ALAN's old one-photo-per-observation canonical cohort. Source-only [run #42446707775](https://github.com/zuizui0223/alan/actions/runs/42446707775) authenticated **597,070** records from one pre-signed original 93.2m photo-index tier, excluded all prior canonical observation UUIDs, and identified **248,071 new observation UUIDs**, including **17,223 on 861 GSHHG islands**. **8,942** new island source observations on **369 islands** already have accepted taxonomy and actual same-year VIIRS/WorldPop source support. Their original photo identities have been frozen into **35 disjoint batches** (34×256 + final 238; [signed batch run](https://github.com/zuizui0223/alan/actions/runs/42447983019)).

The first 256 original photographs on 256 different source-supported islands **have been genuinely acquired and processed** with the unchanged FCP original YOLO11/EfficientSAM ROI-v4: **92 photographs yielded numerical 9-component ROI**; 164 were technically unclassifiable (156 ROI/flip gate, 8 ambiguous palette), and **zero photo downloads failed**. [Original signed image result](results/alan_first256_new_island_real_original_roi9_20261108.json) · [successful source-photo inference](https://github.com/zuizui0223/alan/actions/runs/42447120685). These new photos are distinct observation UUIDs, but one photo/island does not meet the original strict 5-photo/3-species ALAN cell definition. The previous ecological regression remains frozen; the new acquisition is a separately labelled, island-prioritized source cohort, not a globally representative white-flower sample.

## Photo-data expansion: source breadth versus site-year evidence

The SHA-verified original ALAN photo-source queue has **8,942** independent new island observation UUIDs. A newly executed [site-year-focused source selection](docs/ALAN_SITE_DENSIFICATION512_SOURCE_EMPIRICAL_20261108.md) establishes that, after keeping the first 1,280 original-photo targets fixed, the next **512** photos can raise the number of source photo site-years with **≥5 observation opportunities and ≥3 genuine WCVP accepted taxa** from **3 to 60**; the earlier geographic-breadth order would raise them only to **6**. This is an *opportunity*, not already completed colour measurement: **0 new images were inferred by this source-only run**. The source targets **20 islands / 61 actual photo cell-years**, preserves unique photo and observation UUIDs, and remains a separately labeled, nonrepresentative acquisition design. [Completed Actions #42454984788](https://github.com/zuizui0223/alan/actions/runs/42454984788).

## Independent real island photographs — first five acquisition batches

The island-prioritized 93m-source acquisition arm has independently completed original RGB+ROI-v4 inference on **1,280 genuinely new observations**: **492 true numerical original nine-component floral-region ROIs** and **788 technical unknowns**, zero failed original-image downloads, [combined source signature verification](https://github.com/zuizui0223/alan/actions/runs/42456790465). The 1,280 photo/observation UUIDs cover **369 GSHHG islands and 1,060 original photo cell-years**, but just **one** cell has >=5 original photo opportunities and >=3 truly measured WCVP accepted species. This is a different, deliberately island-rich sample from the original global 7,982/710 strict cell-year panel; it is not an increase of 492 regression cells. [Scientific data intake note](docs/ALAN_REAL_FIRST1280_ISLAND_OBSERVATIONS_ROI9_SOURCE_20261108.md) · [aggregate numerical result](results/alan_real1280_original_new_island_ROI9_verified_20261108.json).

## Second 93m source tier: 27,507 further independent island photo opportunities

A new [source-only 1,193,829-photo original tier07 run](https://github.com/zuizui0223/alan/actions/runs/42457757440) has completed SHA-verification of **16 original metadata-only tar groups**. After deduplication against the original canonical and prior tier06 8,942 observation UUIDs, this yields **655,382 further independent observation UUIDs**, of which **50,701** are on GSHHG islands and **27,507** new observations (from **563 islands / 9,493 photo cell-years**) already have accepted WCVP taxa and original-year VIIRS + WorldPop support. Thus the total source-ready, mutually disjoint tier06+tier07 photo observation candidates are **8,942 + 27,507 = 36,449**. **The 27,507 are pre-image-inference source records, not measured white-flower ROIs.** [Source report](docs/ALAN_TIER07_1P2M_REAL_SOURCE_GROWTH_20261108.md) · [verified source counts](results/alan_93m_tier07_1193829_source_only_20261108.json).

## Original tier07 photo ROI9 pilot: 1,024 genuinely new island observations (2026-11-09)

The [first **1,024** original photographs](https://github.com/zuizui0223/alan/actions/runs/42552642250) from the newly authenticated tier07 27,507-observation acquisition source were actually processed using the unchanged frozen FCP YOLO11/EfficientSAM original flower-region nine-palette ROI-v4. Result: **403 genuine numeric nine-component ROI9 successes**, **621 technical unknowns**, zero image retrieval failures, and **57 photo-cell×observation-year units across 29 islands** that genuinely meet ≥5 source photo attempts, ≥3 numeric ROI9 successes and ≥3 distinct successfully photographed WCVP accepted taxa. The source targets span 284 islands; they are *island/site-targeted*, not independent randomized native-flora prevalence data. The original ALAN historical-main five-year regression is unchanged. The initial 70-minute job limit occurred after 896 original photos had been saved; [the successful recovery](https://github.com/zuizui0223/alan/actions/runs/42552642250) measured **only the remaining 128 photos**, without redownloading the first 896. [Source-backed ecological intake note](docs/ALAN_TIER07_FIRST1024_TRUE_ORIGINAL_ROI9_RESULTS_20261109.md) · [verified numeric summary](results/alan_tier07_first1024_actual_originalROI9_20261109_verified.json).

The separately source-disjoint remaining **26,483** original tier07 observation/photo IDs are now being acquired in [52 bounded original-RGB/ROI9 workers](https://github.com/zuizui0223/alan/actions/runs/42553223750) ([Draft PR #13](https://github.com/zuizui0223/alan/pull/13)); future success counts and actual 5-photo/3-measured-taxon sites are **not yet claimed**. The source-level preflight and signed UUID-disjointness checks passed, and no earlier original photo is scheduled again.

## Global 93m-source acquisition: 36,449 real new island photographs / 1,137 strict site-years (2026-11-09)

The **fully measured tier06 (8,942 original photos)** and **tier07 (27,507 original photos)** acquisition cohorts were independently source- and SHA-verified, and [joined in a successful data-only original photo ROI9 audit](https://github.com/zuizui0223/alan/actions/runs/42571444245), with **no image redownload**. The result is **36,449 genuinely distinct original photo/observation UUIDs**, of which **13,842** yield numeric original nine-component flower-region ROI9, spanning **613 source GSHHG islands and 10,724 unique 0.05° site×observation-year units**.

After joining geographically overlapping source sites *without double counting*, **1,137 site×years across 75 islands** actually meet ≥5 source photo attempts, ≥3 numerical nine-palette original ROI9 successes and ≥3 distinct successfully photographed accepted WCVP taxa. The 156 tier06 strict sites and 769 tier07 strict sites overlap at 132 locations/year; **344 additional strict sites became measurable only because genuine photo measurements from both tiers complement each other at exactly the same source island/site/year**. [Source evidence & scientific limits](docs/ALAN_ORIGINAL_ROI9_36449_CROSS_TIER_UNIQUE_STRICT_ISLAND_CELLS_20261109.md) · [inspectable numeric result](results/alan_36449_original_ROI9_unique_strict_island_cells_20261109_verified.json).

**This is additional measurement coverage**, not a causal ALAN response or a native flowering-plant census. Both photo acquisitions were targeted by island and site/year; never append their strict sites unweighted to the previously published historical-main original 710 island site×years. MOIS remains moisture-only.

</details>

湿潤な減少候補の生育地・繁殖形質の類似性と、蛾の花粉運搬が確認された7種について、共通訪花者と同一地点での長期記録率を調べた。[湿地利用のNMDS、GloBI、地点別の比較](results/all_species_forest/seven_moth/README.md)

湿潤な58候補と305比較種を対象に、GloBIの出典をそろえて訪花者とその捕食者の組成を比較した。[全種名簿・生育地をまたぐ相互作用の比較・NMDS](results/all_species_forest/wet_partner_composition/README.md)

ホソヒラタアブとカニグモを含む、湿潤な両群に共通する成虫訪花者・捕食者の組合せを照合し、同地点での記録率と水面距離を比較した。[共通する相互作用候補と水辺の距離効果](results/all_species_forest/wet_shore/README.md)

湿潤な両群について植物の34項目と地点環境の12項目を比較し、別地域での候補種の識別も検証した。[撹乱・遷移・塩分・繁殖・開花期を含む要因探索](results/all_species_forest/wet_drivers/README.md)

淡水湿地と地図で判別できる範囲で、年ごとの水面変化・植生・夜間光と、両群の同じ地点での記録率を照合した。[冠水・水面後退、草刈り後の植生回復と比較可能な地点数](results/all_species_forest/wet_annual/README.md)
