# ALAN三段階メタ解析：PRISMA-2020／PRISMA-S 検索・採択・解析プロトコル草案

**作成：2026-11-19／v0.1／未登録**。過去の肯定・否定の研究結果を確認した後の設計であり、既往の解析を事前登録済みだとは呼ばない。本草案は**今後の検索・採択・抽出・推定条件**を固定するもの。Web of Science／Scopusの検索結果件数は現時点で未取得。

## 研究の中心仮説

ALAN（人工夜間光）の植物繁殖に対する正味効果は、資源・捕食者の再配置と、植物の夜間送粉への機能的依存および昼間・自殖等の代替経路により条件づけられる。この仮説を、異なる論文の3つの効果量を掛け合わせるのではなく、**3段階を別々のアウトカム族としてメタ解析し、同じ照明実験で測定された結合だけを機構の証拠として扱う**。

**段階A（照明→水生起源昆虫・捕食者）**：水生起源昆虫の羽化・流入・捕獲、クモなどの個体数／pitfall活動密度、実際の捕食イベント、人工餌の傷を**別の効果量集合**で解析。主要モデレータは独立に定義した「水域・水辺か、明確に陸域か」。捕獲法・分類群が大きく異なる場合は分け、水辺で必ず正とは仮定しない。

**段階B（照明→訪花・花粉・繁殖）**：植物種別の夜間／昼間訪花率、同種花粉の柱頭への送達、花粉管、果実形成、果実あたり条件付き種子、もともとの花あたり総種子をそれぞれ異なる応答族としてモデル化。主要モデレータは処理以前に決まる独立の植物湿潤性、実際の自主自殖と自家和合性、同じ集団の暗区における夜間**有効送粉**への依存度。自家和合性と自主自殖は同一でない。送粉者の訪花数は同種花粉や種子を直接測った結果の代わりにはならない。

**段階C（補償閾値）**：同じ植物×照明実験の蛾・ハエ等の実測訪花率と、独立研究からの単回訪花あたり同種・和合花粉送達効率qを組み合わせ、植物ごとの「どのqの範囲なら機能的補償が可能か」を条件付きで評価。qが取得できない植物を「補償できた」と分類しない。自殖や花粉制限の有無で種子への変換は異なる。

## 検索対象と実行ログ

探索範囲：各データベースの開始年〜2026-11-19（更新時には更新日を記録）。国・言語・分類群を制限しない。水辺、湿潤、島の語を**必須検索語にしない**（モデレータを持つ研究だけ集める選択バイアスを防ぐ）。Web of Science Core Collection、Scopus、Crossref／OpenAlex、Europe PMC／PubMed、既存メタ解析の文献一覧と引用追跡を別ソースとして保存。Web of ScienceとScopusの検索は**まだ実施していない**。

### Web of Science — Advanced Search

A：照明×水生起源昆虫／捕食者

~~~
TS=(("artificial light at night" OR "artificial lighting at night" OR
ALAN OR "light pollution" OR "night-time light*" OR "nighttime light*" OR
streetlight* OR "street light*" OR skyglow OR "sky glow")
AND (insect* OR invertebrate* OR arthropod* OR spider* OR Araneae OR
predator* OR predat* OR "aquatic emergence" OR "emergent insect*" OR
"aquatic-terrestrial" OR riparian OR "food web" OR "trophic subsid*"))
~~~

B：照明×送粉・植物繁殖

~~~
TS=(("artificial light at night" OR "artificial lighting at night" OR
ALAN OR "light pollution" OR "night-time light*" OR "nighttime light*" OR
streetlight* OR "street light*" OR skyglow OR "sky glow")
AND (pollinat* OR "flower visit*" OR "floral visit*" OR
"pollen deposition" OR "conspecific pollen" OR stigma* OR
fruit* OR seed* OR "reproductive success" OR "reproductive output" OR
moth* OR hawkmoth*))
~~~

C：植物×蛾／ハエの単回訪花効率（**ALANの語は必須でない**）

~~~
TS=(("single visit*" OR "single-visit effectiveness" OR
"single-visit pollen deposition" OR "pollinator effectiveness" OR
"pollen transfer efficiency" OR "per-visit pollen" OR
"conspecific pollen deposition")
AND (moth* OR Lepidoptera OR hawkmoth* OR fly OR flies OR Diptera)
AND (plant* OR flower* OR floral OR pollinat*))
~~~

### Scopus — Advanced Search

A：
~~~
TITLE-ABS-KEY(("artificial light at night" OR "artificial lighting at night"
OR ALAN OR "light pollution" OR streetlight* OR "street light*" OR
skyglow OR "nighttime lighting" OR "night-time lighting")
AND (insect* OR arthropod* OR spider* OR predator* OR predat* OR
"aquatic emergence" OR riparian OR "food web" OR invertebrate*))
~~~

B：
~~~
TITLE-ABS-KEY(("artificial light at night" OR "artificial lighting at night"
OR ALAN OR "light pollution" OR streetlight* OR "street light*" OR
skyglow OR "nighttime lighting" OR "night-time lighting")
AND (pollinat* OR "flower visit*" OR "pollen deposition" OR stigma* OR
fruit* OR seed* OR "reproductive success" OR moth* OR hawkmoth*))
~~~

C：
~~~
TITLE-ABS-KEY(("single visit*" OR "single-visit effectiveness" OR
"pollen deposition" OR "pollinator effectiveness" OR "per-visit pollen")
AND (moth* OR Lepidoptera OR Diptera OR fly OR flies)
AND (plant* OR flower* OR pollinat*))
~~~

**検索品質管理**：ALANが人名などとして引っ掛かる場合も、見出し／要旨スクリーニングで除外。原著で「night lighting」「illumination at night」等を用い検索漏れが分かった場合は拡張式・変更履歴・検索日を追加する。各DBの完全な検索式・検索日・フィルタ・発見件数・RIS/BibTeX原本ダイジェスト・去重前後・二重査読結果・全文除外理由を保存。元レビューの126論文、79論文、PRの9／11候補はPRISMAの「発見件数」ではない。検索未実施はゼロでなく NOT_RUN。

## PECOによる採択／除外

- **Exposure**：人工的な夜間照明（光量・スペクトル・点灯時間・位置）。主要解析は暗対照のある直接照明操作。自然skyglow・衛星光・既存街灯の非無作為化光勾配は別の観察層。人工UV-Bだけの昼夜区別のない曝露は除外。
- **Comparator**：同じ生物・同じ観測努力を有する自然夜光／暗対照。白LED対琥珀灯だけで暗対照がない論文は主要明暗lnRRから除外し、別の色比較階層へ。
- **Organisms/outcomes**：Aでは起源が判明する水生昆虫／捕食者の個体数・活動・捕食イベント。Bでは植物種別に同じ光処理下の夜昼訪花・花粉・果実・種子。模型蛾の色覚のみや人工花行動だけの研究は機能送粉・種子群から除外。
- **Experimental units**：原著の独立灯具設置／物理区画／メソコスムを確認する。各照明群2以上の独立複製を最低条件とし、種・個体・花・夜間トラップ反復を灯具複製に増やさない。Lizard30の植物・クモ2論文、Penryn54の栄養段階3論文、同一可能性のあるLuの2論文は独立研究として加算禁止。
- **Extractability**：効果量と分散、独立ユニット数、対照と処理、観測努力、変動の測定・推定が必要。分散不足は資料請求候補／記述プール、p値のみから偽の標準誤差は生成しない。「有意でない」は効果量0でない。
- **Moderators**：水辺／陸地は照明導入前の生態系構造、植物湿潤性は独立した生育地属性（地域に妥当なEllenberg等）、自家和合性・自主自殖・送粉効率は原研究か独立外部の事前データから確認。欠損は unknown とし、乾燥・非自殖・低効率に勝手に代入しない。
- **Studies**：原著、プレプリント（感度層）、未公表博士論文／技術報告も探索して出版バイアスを弱める。レビュー・同一実験の二重出版は引用発見用であり独立元研究に加えない。
- **Risk of bias**：無作為割付、明暗区の初期条件、真の設置数、ブラインド評価、観測努力、捕獲法（pitfall活動密度≠個体数）、UV波長、昼夜のズレ、夜間点灯と日没の交絡、反復の扱い、結果の選択的報告を独立採点。観察・操作群は主要解析を混合しない。

## 効果量と多層モデル

**符号規則**は常に「照明区／暗区」。正の連続平均・観測努力が等しい条件で lnRR = log(mean_L / mean_D)、独立照明区画の近似分散は SD_L²/(n_L mean_L²) + SD_D²/(n_D mean_D²)。nは独立灯具設置／区画の数で、捕獲昆虫や花の総数ではない。

**二項果実成功**は実際に曝露した花を分母とする log risk ratio、Poisson計数は独立調査努力をoffsetとした log rate ratioを別familyで扱う。攻撃痕のlog odds ratioは log RRと同じでない。0、非正値、過分散、共有対照、対応する照明ブロックでは原著設計に基づくGLMMと共分散が必要。恣意的に0.5を足すのは事前指定した感度分析に限る。

Sanders 2020／2021元データの1,304効果量は**Hedges d**、Guenat & Dallimer 2023は**Hedges g**であり、原腕平均・分散・標本数を回復せず lnRR に読み替えない。

各エンドポイント別に metafor::rma.mv を使い、真の source programme／paper／physical experiment／effect familyを識別する。Vには共有対照や同一照明区・同一種の共分散を可能な限り保存し、未知共分散は異なる仮定の感度解析。study／programme単位のクラスタ頑健推定を併用。原研究の p値・効果量方向のみをメタ分散に作り変えない。

**Aの一次モデレータ：水辺 vs 陸地**。**Bの一次モデレータ：独立した湿潤生育地嗜好**。加えて光色・照度・夜間点灯率・捕獲法・自主自殖・暗区の実測夜間有効送粉依存は事前指定するが、研究数が不足したら単変量・記述に降格。15〜20研究は計画開始の目安であり、モデレータ検出力の根拠ではない。各カテゴリ数と外挿可能性・study間の共通支持を検証し、シミュレーションで必要なkを評価する。

研究数が少ないときのファンネル／Eggerは自動実施しない。元データの選択・アウトカムの報告バイアスや小研究効果は、十分な独立programmeと比較可能な推定対象があるときだけ階層構造に対応した方法で探索する。

## 機能的補償（C）の識別条件

植物iの有効花粉受領の代理 Φ_it = Σ_g V_igt q_igt とする。Vは植物別の花あたり／時間あたり訪花率、qは同植物の適合する柱頭への**和合花粉**の単回訪花効率。蛾の訪花がΔV_moth<0、ハエがΔV_fly>0で、qが照明に依存せず、他ギルドが不変なら、必要な効率比は q_fly/q_moth ≥ −ΔV_moth/ΔV_fly。分母がゼロまたは負ならそのギルドによる補償は成り立たない。自殖や花粉制限の非線形を別途確認。qが取得できない植物を補償「成功例」に数えない。「補償できる植物割合」は独立qと同じ植物のVと不確実性が揃ってからであり、現段階は必要効率比・部分識別範囲のみ。

## 独立した湿潤植物マクロ主結果との接続

mainの長期減少候補の湿潤性偏り β≈+0.534 は**植物種群の構成差**。ALANによる植物の因果的死亡・送粉喪失の推定値ではない。別の英国NPMS実測被度の湿潤性×光の係数は強い季節・生育地・調査区分感度でゼロを跨ぐ結果がある。A/Bの一次モデレータが独立実験で再現されなければ、主結果の「機構的証明」を宣言しない。逆符号・無解像も明記する。

## 既存研究と作業開始時の母集団

機械可読seed台帳：config/alan_prisma_three_stage_seed_census_20261119.json。Aに**11件の意図的候補**（作成時の9件＋2025水生双翅目・2017同一水辺プログラム）、Bに**11件の意図的候補**（観察や訪花だけ、重複可能性含む）。これは検索結果の網羅性・採択数ではない。Seymoureの79論文はreviewの独立原著、Guenat照明22候補（UV-B8、非UV-B14候補）は原著重複を含み、Sandersの126論文／1,304効果量は元のHedges dデータセットであって、これらを単純加算しない。

既存監査の**12植物繁殖原著（9直接光照射・3観察）**で、独立した暗所における夜間送粉依存と、同じ植物の明暗全日開放種子効果＋SEが完全にそろった研究は0。**これはBの一般的な訪花lnRR研究が0という主張ではない**。すべて再抽出が必要。

Sanders Dryad DOI 10.5061/dryad.wpzgmsbjnの公開CSV（Sanders_et_al_datafile.csv）と解析RコードはZenodo原ファイルのMD5/SHA256により真正取得済み。著者のKeep 1,109／Delete 195を分離し、AではKeep152効果・23元Study（うち野外Keep130・19Study）、Bの植物繁殖関連ではKeep12効果・1Studyと再集計した。段階Aの記録は水辺・陸地／羽化・流下・捕獲を同じエンドポイントに自動昇格できず、[実測監査の追加結果](ALAN_THREE_STAGE_SANDERS_KEEP_DELETE_AQUATIC_EFFECT_STAGE_FULLTEXT_TRIAGE_20261119.md)に保全した。次にGuenat2023真正S4/S7/S13、Seymoure79件の原著DOIを重複照合し、2021〜2026の新照明実験を追加する。表層の信号変化だけでなく、原実験の暗対照・n・SE・水辺・湿潤・自殖の原情報を集め、PRISMA flow実際の検索日・nを別ファイルに保存する。

参照：Sanders 2021 DOI 10.1038/s41559-020-01322-x（元Dryad 10.5061/dryad.wpzgmsbjn）、Guenat & Dallimer 2023 DOI 10.1002/ece3.9990、Seymoure 2023 DOI 10.1098/rstb.2022.0356、PRISMA 2020 DOI 10.1136/bmj.n71、PRISMA-S DOI 10.1186/s13643-020-01542-z、metafor公式文書。


## 2026-11-19 真正元データからの採択品質ゲート

[Sanders2021元作者アーカイブ実行CI](https://github.com/zuizui0223/alan/actions/runs/44529932726) では元CSV1,304行・元著者RコードのMD5/SHAを確認。さらに[元著者Keep/Deleteを段階別に分割したCI](https://github.com/zuizui0223/alan/actions/runs/44538857424) は4テスト成功、[段階Aの物理サイトと測定過程を分離するCI](https://github.com/zuizui0223/alan/actions/runs/44540344752) も成功。公開CSVの「Study=23」や「n=153候補」は**適格な独立照明研究数ではない**。

水際のManfrin2017と2025は同一Westhavellandの照明1地点・対照1地点に入れ子であり、追加引用を独立した12灯×2研究と数えない。4森林河川ペアのPerkin2014では水生無脊椎動物のドリフトは減ったが羽化は不明瞭、Westhavelland2014では水生Dipteraの羽化量が減り灯下飛翔捕獲は増えた。**捕獲法・移動・羽化の区別は水辺モデレータより先に必要**。これらは事前登録ではなく、検索前の原著検証として扱う。
