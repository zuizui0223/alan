# ALAN: 島嶼の「夜間送粉への進化的シフト × 長期光害」— 同地点ペアの実測検出可能性（2026-11-10）

## 判定

**厳密な3重交互作用を推定できる独立植物ペアは、現有の写真・独立送粉資料では確認できなかった。** したがって三重交互作用の β、信頼区間、p値は算出しない。「島で夜間送粉に乗り換えた」ことは、単に島の植物が夜咲き・淡色・GloBI蛾陽性であることでは証明できない。祖先・姉妹系統の有効な昼行性送粉と島系統の有効な蛾送粉を独立資料で裏付ける必要がある。

### 真のサンプル取得元

- ALAN [元93,209,351写真ID、GSHHG正確な島、WCVP accepted植物](https://github.com/zuizui0223/alan/actions/runs/42616279235)、16独立の原本source group。各植物種×実座標0.05°site×yearのsource `photo_ids` と独立`observation_uuids`を積み上げた（写真同士は生物学的独立ではない）。
- 2019・2021・2023・2024・2025年の記録地点を、[元33年PCNL衛星夜間光](https://github.com/zuizui0223/alan/actions/runs/42624285048)（1992〜撮影前年まで）と結合。元85,318地点×年、島内・同年に光量差がある全562 GSHHG島の77,647地点×年に絞った数も**別カラム**で計算。
- [Kendall & Nicholson 2025、昼夜送粉除外実験](https://doi.org/10.1111/ele.70036)の元データ・種形質・研究報告のraw Git blob SHAを固定した[PR #29](https://github.com/zuizui0223/alan/pull/29)。published最終139種と生データの145種名は同一と仮定しない。過去データで**54種**がWCVP accepted実測写真に正確に対応。
- 上記54種に、少数の**文献確認済み島固有または対照候補（seed）**を付加。これは**全球の島固有夜間送粉植物の網羅的な文献目録ではない**。

## 1. 訪花／生物学的な送粉を直接検証した島植物の源写真サポート

| 植物と証拠の強さ | 元写真ID | 独立観察UUID | GSHHG島単位 | 原本source地点×年 | 562島PCNL解析に含まれる地点×年 | その種が同島・同年の複数明暗地点で記録された島 |
|---|---:|---:|---:|---:|---:|---:|
| **Ipomoea habeliana** — ガラパゴス固有、**野外でスズメガが唯一の有効送粉昆虫**と記録。自家受粉も可能。日→夜への祖先状態の進化的移行は不明 | **44** | **22** | **5** | **10** | **3** | **0** |
| **Brighamia insignis** — ハワイ固有、蛾媒型の花形質。直接有効蛾送粉は未確認、野生絶滅（IUCN）。進化的移行不明 | 0 | 0 | 0 | 0 | 0 | 0 |
| **Canarina canariensis** — カナリア固有、鳥の実効送粉が実験検証された対照候補 | 917 | 472 | 5 | 128 | 126 | 4 |
| **Digitalis canariensis**（旧 Isoplexis canariensis）— カナリア固有、鳥の実効送粉。**WCVP accepted synonym解決済** | 396 | 194 | 2 | 55 | 54 | 1 |
| **Echium simplex** — テネリフェ固有、飛翔昆虫の送粉貢献が確認。**昼限定ではない** | 121 | 56 | 4 | 25 | 25 | 2 |
| Cordia lutea — ガラパゴス在来だが**群島固有ではなく**、昼夜両方の送粉／訪花証拠。昼行性だけの対照には**使えない** | 1,000 | 605 | 19 | 114 | 27 | 3 |
| Ipomoea carnea — 広域分布、実験研究の記録では昼夜開花の候補。I. habelianaとは生活形も違い、同一セル×年での共記録は0 | 347 | 160 | 46 | 134 | 93 | 7 |

原本sourceには、異なる島ポリゴンの栽培・移植観察を区別する`captive/cultivated`フラグがなく、すべてのGSHHG島の写真がその種の*自生*を証明しない。例えば `Echium pininana` は31 GSHHG島に元写真があるが、31島に自生する島嶼固有種という意味ではない。**単に島で見つかった分類群を「島固有」と名付けてはいけない。**

## 2. 実験証拠に基づく「同地点ペア」の全数監査

Kendall & Nicholson の元昼夜送粉除外実験の種で、ALANの562島に元写真記録がある**54種**の間で：
- 同一GSHHG島×0.05°セル×観測年に両種が記録された**298組合せ**。
- そのうち植物属が同じ組合せは**8組**（Silene、Lilium、Asclepias）。
- その8組のうち、原著形質の`bloom_period_simple`が一方 **night**・もう一方 **day** と判明する**同属ペアは0組**。夜咲きと昼咲きは**有効送粉への依存**の強さを意味しないので、これはさらに弱い仮の前段階のスクリーニングである。
- *Ipomoea habeliana* と *Cordia lutea* は**3島の8セル×年**で共記録がある。しかし *C. lutea* には夜間の蛾による花粉運搬も報告され、属・生活形・色も一致しない。**昼行性送粉のみの負の対照には不適格**。
- *Ipomoea habeliana* と実験形質情報のある *Ipomoea carnea*（昼夜開花、低木）は、共通GSHHG島は1つでも、**同一セル×年の観察が0**。*I. carnea* は同属・昼限定の適合対照ではない。
- *Digitalis canariensis*と*Canarina canariensis*はカナリアで比較できる鳥送粉の**対照候補同士**であり、ガラパゴス夜間蛾送粉との同じ島の負の対照にはなれない。

## 3. 引き続き検証すべき仮説と順序

元の「ALAN × 夜間送粉への進化的シフト × 島／本土」は魅力的だが、現時点で直接推定するのは生態学的に誤認を招く：
1. **独立系統の送粉者・祖先状態証拠**を揃える。カタログの「蛾訪花陽性」ではなく、雌しべへ有効花粉を届ける種を分け、昼の代替・自殖率も定量。
2. **島ごとの自然植生か／植栽かを独立に検証**。WCVPのTDWG3地域在来区分は、各島自生とは別もの。検出可能性と花の撮影魅力度の差は、同一地点比だけで必ずしも相殺されない。
3. **自生地点の複数明暗サイトと、真の昼行性・夜行性比較種**を揃える。同じ島・季節・観察努力で両種の観測機会が存在することを必要条件にする。
4. 島／本土を追加するなら、**同じ島内sourceだけでなく本土の比較母集団を独立整備**し、栽培・外来混入の程度と種ごとの検出確率が共有可能であることを検証。原本93mの中の本土標本を無条件で島標本の対照にしない。

**結果を都合よく正・負の回帰係数として生成しない。条件未達なら検定不能として保存する。** 本PRは実際の元写真IDと独立送粉文献に基づく**検出力／比較サポートの診断**であり、「光害影響なし」という結論ではない。

## 外部出典
- [McMullen 2009: Ipomoea habeliana 夜間スズメガ送粉](https://doi.org/10.1111/j.1095-8339.2009.00963.x)。
- [Brighamia insignis: 蛾媒型だが実効蛾送粉は不明](https://doi.org/10.1016/j.flora.2019.151448)。[IUCN 野生絶滅](https://doi.org/10.2305/IUCN.UK.2023-1.RLTS.T44080A184207893.en)。
- [Kendall & Nicholson 2025 の昼夜送粉除外メタ解析](https://doi.org/10.1111/ele.70036)。
- [Canarina を含むカナリアの鳥送粉](https://pubmed.ncbi.nlm.nih.gov/18931830/)、[旧 Isoplexis（現 Digitalis canariensis）鳥送粉](https://pubmed.ncbi.nlm.nih.gov/21628148/)、[Digitalis canariensis のaccepted名](https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:801983-1)。
- [Echium simplex の送粉実験](https://pubmed.ncbi.nlm.nih.gov/30334319/)。
- [Cordia lutea と夜間花粉運搬・ガラパゴス在来](https://pmc.ncbi.nlm.nih.gov/articles/PMC3489146/)。

## 再現
[元写真・昼夜送粉実験の独立証拠×同一セルペア全数監査スクリプト](../scripts/analysis/alan_island_endemic_night_moth_pair_support_20261110.py)、[元原本のsource-gated Actions](../.github/workflows/alan-island-nightshift-pair-source-feasibility.yml)、[元データ再現ジョブ](https://github.com/zuizui0223/alan/actions/runs/42747528696)、[結果の正式アーティファクト](https://github.com/zuizui0223/alan/actions/runs/42747528696/artifacts/13732565503)。
