# ALANの肯定結果・反証・生物間相互作用の統合（2026-11-18）

## 統合したい生態学的発見

**人工夜間光（ALAN）の生態学的影響は、動物が増えた／減ったという総数だけでは決まらない。照明の周囲では水域由来昆虫と捕食者が集積し、花では夜行性チョウ目からハエ目への訪花経路の再配分が生じうる。植物が実際に利用できる花粉送達経路は、見かけの訪花ネットワークより少ない傾向があるため、別の訪花者の増加は繁殖機能の代替を保証しない。**

これは互いに独立した原実験と観察ネットワークを**段階別に統合した解釈**であり、湿潤植物について ALAN→水生昆虫→捕食者→送粉→結実 の全媒介経路を因果的に同定したという主張ではない。

元mainで固定されている唯一の主要マクロ結果は、増光地点での長期記録率減少候補に**湿潤な生育地を好む植物が過剰**という種群構成である。別々の生態系で得られた実験結果を、この種群偏りの「確定機構」に昇格させない。

## 証拠の階層（独立単位を固定した原本の結果）

| 段階／証拠の位置づけ | 観測単位・元データ | 肯定結果 | 反証・限界 |
|---|---|---|---|
| A. 元mainの種群構成（観察マクロ） | 1992–2024の3,901評価種から283減少候補；環境調整は158候補と1,419参照植物 | 湿潤性差 **+0.534**、主要調査地域ジャックナイフ近似95%区間 **[+0.158,+0.910]**。初期の湿潤性平均差 **+0.59**、17指標BH **q=.030** | 12地域の7地域で正・5地域で負。記録率は存在量・局所絶滅・結実ではない。減少候補は事後スクリーニング |
| B. 水域起源の餌資源→捕食者（ALAN野外操作） | Parkinson et al. 2020、**明10＋暗10＝20独立メソコスム**、クモの3時点は反復測定 | 水生起源昆虫のpan捕獲 **6.63倍 [4.29,10.11]**、アシナガグモ類 **2.01倍 [1.66,2.44]**、網獲物 **2.40倍 [1.67,3.67]**（いずれも20区画bootstrap 95%区間） | 陸生起源昆虫は **1.30倍 [0.80,2.03]**。網獲物／クモは **1.19倍 [0.81,1.82]**。蛾・花粉・結実効果は実測なし |
| C. 夜間訪花者ギルドの再配置（ALAN野外操作） | Lu et al. 2026、**5明暗ペア**、元表400集計／1,038焦点植物／1,606ネットワーク（有効1,605行、44植物ラベル） | 同一の群集集計表で夜間チョウ目 **暗838→明357**、夜間ハエ目 **暗54→明151**。5ペアとも夜間チョウ目／ハエの相対比低下 | 5ペアの厳密両側符号反転 **p=.0625**。3表の訪花分母は一致せず、統合補償率を作らない |
| D. 同じ植物上の時間帯・ギルド交代（ALAN野外操作内） | Lu原ネットワーク、昼夜とも比較できる **33植物×ブロック組（14植物、5独立照明ペア）** | 夜間チョウ目↓の21組で夜ハエ↑11、昼ハエ↑12、どちらかハエ↑**15**、両方↑8 | どちらの時間帯でもハエ↑なし**6**。再訪問回数を独立処理区として数えず、和合花粉量・種子機能は未測定 |
| E. 機能的な繁殖の違い（別のALAN実験論文） | Lu et al. Journal of Plant Ecology 2026、焦点植物の繁殖・受粉 | 自家不和合 *Pedicularis kansuensis* の夜間同種花粉受領 **−53.8%**、種子成績 **−56%**。自家和合 *Euphrasia regelii* では明瞭な種子低下は検出されず | Luの別論文と研究地・著者が重なる。個体・照明区の重複を確かめるまでは独立なメタ研究2件と数えない。湿地への直接外挿不可 |
| F. 訪花頻度≠花粉送達（非ALAN観察マクロ） | Page et al. 2021 Figshare原623行から**33研究・72植物×研究×期間ネットワーク／62学名** | 訪花割合vと同種花粉SVE×訪花率割合pの全変動距離 **平均0.2098 [0.1604,0.2638]**。機能的Hill2−頻度Hill2 **−0.465 [−0.615,−0.240]**、**50/72**で機能的経路が狭い | 研究IDクラスタの記述的区間であり、照明操作の効果量・植物適応度・送粉者喪失後の損失量ではない |
| G. 同じ植物の蛾とハエの機能比（非ALAN観察） | Page元CSVの同指標・同単位で比較できる**5研究5植物** | 相対花粉送達率 q×V（ハエ／蛾）は *Amomum subulatum* **1.009**、*Erythranthe lutea* **0.00720**、*Silene caroliniana* **0.0765**。単回効率と訪花率の積が重要 | 別の1植物は訪花率単位が不一致、もう1植物は両者q=0で比不能。ALAN実験でも湿潤種でもなく、蛾・ハエの全球平均効率を推定しない |
| H. 希少訪花者の機能的寄与（非ALAN観察） | Ballantyne et al. 2015、英国**1ヒースランド群集**、5植物／1338原訪花 | 湿潤減少候補の *Erica tetralix* で *Bombus jonellus* は **11訪花（2.76%）** から **608同種花粉（12.99%）** を送達 | 当該マルハナバチは4観測日。5植物とも最多訪花者と最多花粉送達者は一致。蛾訪花0、ALAN操作なし |
| I. 湿潤植物の訪花者→捕食者（文献記録の鎖） | 元の**363植物レコード／361ユニーク学名**；記名訪花者72・捕食者143、272文書上の鎖 | 捕食者の可能なギルドはクモ・鳥・トンボ・捕食昆虫・カマキリ・爬虫類・両生類に及ぶ | 272鎖の両端でiNaturalistを含まない出典組は**2のみ**。水辺植物―蛾―クモの8件も捕食側すべてiNaturalistで、同所の捕食／ALAN／送粉を示さない |

**統計上の単位**：20メソコスムと5明暗ペアを「60クモ測定」「1,605訪花データ」として擬似反復しない。33元研究内の72ネットワークも72独立研究とは呼ばない。元のParkinson 4反応は同一実験の別endpointで、4件の独立ALAN実験ではない。

## 一つの理論に接続する：夜間光が「利用可能な相互作用経路」を組み替える

1. **資源の再配分**：水域由来の昆虫が照明周辺の捕獲地点に集中し、クモ個体数と網獲物が同時に増える。ただしこの資料から捕食による送粉者除去効果は計算できない。
2. **訪花の再配分**：夜間チョウ目への依存が減る代わり、ハエなど別ギルド・時間帯への訪花が一部植物に生じる。同じ植物・同じ場所で代替訪花がない場合もある。
3. **機能的なボトルネック**：見かけ上の訪花回数 (V_{ig}) と、単回の同種・和合花粉送達 (q_{ig}) は別の量。植物 i が受け取る花粉の近似は (P_i ≈ \sum_g V_{ig}q_{ig})。ある訪花者が増えてもqが小さければ元の経路は補えない。
4. **繁殖上の保険**：自動自殖や昼間の実効送粉、脊椎動物の送粉などが実際に有効なら、夜間依存に由来する繁殖損失を弱める可能性がある。これは全植物に確認された実証ではない。

同じ植物・照明ブロックで夜間蛾の訪花差ΔV_moth<0、ハエ訪花差ΔV_fly>0のとき、他経路の変動をゼロとする**仮想的な必要条件**は q_fly/q_moth ≥ −ΔV_moth/ΔV_fly。Lu原本の植物例は *P. kansuensis* **0.67以上**、*Anaphalis hancockii* **23.5以上**、*Carum carvi* **2.75以上**。qは未観測で、これは補償成功率ではない。

**最大の非同定箇所**：Parkinson(水辺・餌・クモ)、Lu(高山草地・訪花)、Page/ Ballantyne(非照明・機能花粉)に同じ植物×同じ照明区×同じ捕食者／送粉者のデータがない。ALANの湿潤植物減少を「捕食者が媒介した」「ハエが受粉を救った」と結論しない。

## 競合説明と不採用の「肯定的」解釈も保持

- 元mainの湿潤性調整差は**種の選抜・観測分布の記述**。調査者、地域、対象種、湿地の水位・排水・植生管理の交絡が残る。
- [PR #86](https://github.com/zuizui0223/alan/pull/86) のNPMS追加感度では、湿地×種湿潤性×夜間光の交互作用が調査record_statusを一定にした系列で**−0.861 [−1.778,+0.057]**、生育地安定＋記録区分安定の共同制限では**+0.200 [−1.125,+1.525]**（二方向CI）となる。したがってこの三重交互作用を新たな主要「肯定結果」としてマージしない。
- 元mainの訪花者・捕食者組成比較では、候補植物間の共有相手種が特有に多いという6検定の**最小BH q=.272**。記録の鎖を同所実験の因果ネットワークと誤認しない。
- [PR #85](https://github.com/zuizui0223/alan/pull/85) の夜間繁殖依存群×島の元写真記録率の厳格22種比較は**−0.368 pp [−2.025,+1.302]**でゼロを跨ぐ。別の [PR #84](https://github.com/zuizui0223/alan/pull/84) でも自動自殖のALAN緩衝効果は独立島に対して頑健ではない。写真記録率を結実や送粉機能に置換しない。
- Pageの多ギルド55ネットワークの平均TVは単一ギルド17ネットワークより大きいが、元研究内に両群あるのは**1研究のみ**。ギルドの多さが原因で機能の乖離が増えるという係数は**推定停止**。
- PageとBallantyneには**同一の元研究**が含まれる。Földesiの原メタにもBallantyne由来の効果量が含まれ、独立観測として二重計上しない。
- Parkinsonの原資料にはpan捕獲総数1992と論文の1990、クモ数の原集計比と論文表現に不整合が残る。元値の修正・出典単位の混合は行わない。

## 研究の暫定主張（原稿で言えること）

> Artificial light at night is associated with a moisture-affinity bias among plants showing declining long-term recording rates. Independent lighting experiments demonstrate redistribution of aquatic-origin insect captures, spider abundance, and nocturnal versus diurnal flower visitors. A cross-study functional pollination synthesis shows that visitation networks overstate the number of effective pollen-delivery pathways. Whether these interaction changes explain the moisture-affinity bias or reduce seed production remains unresolved.

**主結果は湿潤性の種群偏り、機構支持は複数独立研究の段階的な証拠、統一媒介効果は未同定**とする。この区別が、肯定結果を失わず因果の言い過ぎを避ける最も安全な統合。

## 出典・再現可能性

- 主結果：[元main原稿](../manuscript/ALAN_WET_HABITAT_ja.md)、[湿潤性環境調整](../results/all_species_forest/wet_environment/README.md)、[湿潤種相互作用の否定も含む探索](../results/all_species_forest/wet_partner_composition/README.md)。
- 原本20メソコスム：[統計・出典監査](ALAN_BIOTIC_INTERACTION_SOURCE_20261118.md)、[固定JSON](../results/all_species_forest/biotic_interactions/parkinson2020_source_verified_20261118.json)、[原実験再解析成功CI](https://github.com/zuizui0223/alan/actions/runs/44215675461)。
- Lu 5明暗ペア：[同花・昼夜の訪花者交代](ALAN_BIOTIC_LU2026_DAYNIGHT_FUNCTIONAL_SUBSTITUTION_20261118.md)、[固定JSON](../results/all_species_forest/biotic_interactions/lu2026_fiveblock_sameflower_temporal_rewiring_verified_20261118.json)、[原本成功CI](https://github.com/zuizui0223/alan/actions/runs/44245596409)。
- Page 33研究72ネットワーク：[機能的ボトルネック](ALAN_BIOTIC_PAGE2021_FUNCTIONAL_BOTTLENECKS_20261118.md)、[同一植物での蛾・ハエの真のq×V](ALAN_BIOTIC_PAGE2021_MOTH_FLY_SVE_FLUX_20261118.md)、[固定JSON](../results/all_species_forest/biotic_interactions/page2021_72network_functional_bottleneck_source_20261118.json)、[真正原本成功CI](https://github.com/zuizui0223/alan/actions/runs/44268234628)。
- Ballantyne / Földesi：[機能ネットワークと重複監査](ALAN_BIOTIC_FUNCTIONAL_POLLEN_REDUNDANCY_HEATH_20261118.md)。
- 研究間独立性：[14行のendpoint別出典台帳](../results/all_species_forest/biotic_interactions/experimental_stage_registry_20261118.csv)。**14行は14独立研究ではない**。

この統合は「元mainの長期湿潤性結果を書き換えない」「暫定的な相互作用候補の根拠を検証可能な原データ・コードとともに保全する」という目的でPR #87からmainへ移す。