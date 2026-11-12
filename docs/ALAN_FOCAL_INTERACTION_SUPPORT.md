# 大きな正負の寄与を示した6種の相互作用と季節の観測範囲

共通544種から島−本土の符号付き寄与の上位・下位各3種を選んだ。元の観測年は2019、2021、2023、2024、2025年、全体は78地域・195地域年・10,698地点年である。植物への夜間光係数は再推定せず、同じ種でどの機構が現地記録に結び付くかを調べた。[6種の機械可読表](../results/focal_species/interaction_support.csv)・[抽出コード](../scripts/analysis/alan_focal_interaction_support.py)

| 植物 | 名前でつながる捕食者 → 成虫の蛾 | 3者が同年の1°にそろう元の地点年（島／本土） | 開花状態調査（島／本土、件） | 努力量付き訪花調査（島／本土、分） |
|---|---|---:|---:|---:|
| *Asclepias syriaca* | 該当記録なし | 0／0 | 70／6,047 | 0／1 |
| *Trifolium repens* | *Argiope bruennichi* → *Euclidia glyphica* | 0／6 | 4／377 | 290／56 |
| *Stellaria media* | 該当記録なし | 0／0 | 0／0 | 0／10 |
| *Cirsium vulgare* | 該当記録なし | 0／0 | 0／110 | 20／0 |
| *Hedera helix* | *Hemidactylus frenatus* → *Spoladea recurvalis* | 0／0 | 0／0 | 0／0 |
| *Dittrichia viscosa* | *Hemidactylus frenatus* → *Spoladea recurvalis* | 0／0 | 0／0 | 0／0 |

**捕食の地理的支持。** 既存のヤモリ2種、カマキリ、コガネグモ科の *A. bruennichi*、ツマアカスズメバチの記録について、捕食資料での成虫同定、蛾の種名、旧GloBIにおける同じ蛾の訪花植物の種名を要求した。GBIFの陽性人為観察を年・物理陸塊・格子で照合し、蛾の現地記録にも成虫ステージと座標誤差1 km以下を要求した。大きな格子でそろうシロツメクサの本土6地点年も、元の0.05°でそろう地点年は0。6種すべてが細かい格子で島0／本土0だった。各数値は条件を満たした元の植物地点年の数であり、異なる日の動物・植物観察も含む。捕食回数、同じ花の前での攻撃、被食率は測っていない。0はこの資料の一致記録がないことを表す。[昆虫・クモの成虫捕食出典と地理集計](../results/alan_arthropod_adult_moth_flora_20261111.json)・[ヤモリ捕食記録](../results/alan_gecko_primary_fiveyear_20261111.json)

**捕食者が送粉も担う可能性。** 成虫の蛾を捕食し送粉もすると分類された80種のうち島記録をもつ捕食者では、一次記録の再確認後、捕食者と蛾が同じ植物に結び付く真正の捕食候補に *Heracleum maximum* が残った。焦点6種にこの条件を満たす例は0種だった。例えばキヅタへ花粉を運ぶスズメバチがいても、キヅタを訪花する蛾をそのハチが同じ場所で捕り、失われた送粉を補ったことはまだ確認できない。[一次記録まで確認した種の組合せ](../results/dual_role_80/island_cohort_reviewed_links.csv)・[コホートと空間支持](../results/dual_role_80/summary.json)

**訪花調査の分母。** 表はUK PoMSの10分定時調査と、Budburst／UGAの元の植物・地域年に一致し、正の観察時間と解釈可能な非鱗翅類の広域区分を持つ調査を対象にした。無訪花と記録された観察も時間に含めた。英国は種の二名法に直接一致する標準化対象名を用いた。陸塊が判別できた調査のみを集計した。シロツメクサの島290分は英国36_28の2021年120分・2023年170分、本土56分は米国16_23・20_25の2023年46分・2024年10分。6種とも同じ地域・年の両側で同じ種を比較できた例は0。トウワタの2024年の本土にはほかに2分調査があるが分類群別の内訳が不明なため、表に含まれない。ドイツのEuPPollNetにもシロツメクサとコハコベの訪花記録があり、比較可能な植物別観察時間を確認できなかった。[英国原データ](https://doi.org/10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618)・[調査源別の植物照合](../results/shared_visitors_regional/regional_plant_support.csv)・[北米のイベント原表](../results/bilateral_visitors/us_candidate_events.csv)

**フェノロジー。** 開花状態の件数はUSA-NPN、スウェーデン、豪州ClimateWatchの対象地域年における開花・未開花の明示記録である。元の写真地点から離れた同地域年のモニタリング地点を含む。*A. syriaca* は島で開花7/70、本土433/6,047の記録があり14地域年にまたがるが、同じ地域年の島・本土×明・暗4区分のうち最多でも3区分。*T. repens* は島0/4、本土96/377、最大2区分。*C. vulgare* は本土4/110、最大1区分だった。件数の比率は観察季節・頻度に依存する。どの種にも同地域年の4区分に開花記録がそろわず、局地的な蛾の活動曲線と花期の重なりの変化量は未推定である。[モニタリングの抽出範囲と四区分](../results/regional_repeated_phenology/README.md)

再計算には、既存のデータ取得時に保存した作業フォルダを渡す。

```
python scripts/analysis/alan_focal_interaction_support.py WORK_DIR results/focal_species/interaction_support.csv
```
