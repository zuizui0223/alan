# 湿潤植物偏りの機構的Resultsを、独立ALAN操作実験の植物形質で再現できるか

2026-11-20 | Anic (Penryn-54) + Giavi (Swiss-6) original-source audit
**Exploratory effect modification only — no wet-plant causal mechanism identified.**

## 主結果は変更しない

alan/main の植物3,901評価種の大規模記録率解析では、増光に関連して記録減少候補に選ばれた植物は比較植物より湿潤性が高い（調整後+0.534、地域区間は正）。この観察的な偏りは維持するが、同一種×同一地点×同一年の ALAN 因果効果を示したわけではない。

本資料の新しい問い：「別の実験で、植物種の湿潤性が高いほど、照明操作に対して植物資源量や送粉者訪花が悪化するか？」 二つの**独立した元照明設備**があり、元実験を混ぜて合成した単一のメタ係数はない。

## A. Anic et al. 2022 — 植物バイオマス

Anic, Gaston, Davies, Bennie (2022, DOI 10.1002/ece3.9157) の原 Figshare PDF（SHA256 baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0）から、**13植物種、665植物×元トラフ行、54独立物理メソコスム、9セル×各6トラフ**の原平均・SDを再構築している。

英国Kewの公的 BIFloraExplorer source BI_MAIN.csv（Git blob SHA1 e6a1da929af81f370d2a1b79177b489e3e573a06；英国 PLANTATT Ellenberg moisture F）を原種学名で照合。**12植物は二名法完全一致、Agrostis tenuis はKew POWOで正式に Agrostis capillaris の同物異名であるため1例のみ明示的に正規化**。元13種すべてでFが得られ、F=4:2種、F=5:8種、F=6:2種、F=8:1種。

植物ごとに元3栄養条件の等重み平均総乾燥バイオマスを取り、L2/L1, L3/L1, L3/L2 の3つの光処理コード間のlog比を比較。3対比共通完全例12種を固定しFと効果のSpearmanの絶対最大値を事前に比較した4999回の種ラベル置換帰無で多重比較補正。**重要：原 Figshare のL1/L2/L3コードと暗・琥珀・白色LEDの著者キー対応は未照合。原著の「白色LEDでLotus pedunculatus減少」を符号マッピングの根拠に逆用しない。**

| 元の光コード比 | 同じ12種Spearman rho | 無補正種ラベル置換 p | 3比較最大統計量補正後 p |
|---|---:|---:|---:|
| L2/L1 | **−0.6344** | .0308 | **.0766** |
| L3/L1 | −0.5052 | .1040 | .2358 |
| L3/L2 | −0.0352 | .9276 | .9992 |

L2/L1の負方向は最大F種Lotus pedunculatusを除いてもrho−0.5721で、単一種だけが符号を決めているわけではない。しかし、F>7の植物は13種中1種のため、湿潤種一般対乾性種の交互作用とは解釈できない。また13植物の生産量には同じ54元トラフの共有共分散があり、植物種を13独立照明実験と扱ったp値ではない。元Anic照明コードと対照がどの組か未特定、元2018植食者・2015花数との54トラフIDは未結合。育成用の土壌水分操作もない。

[真正元データ＋Kew F原文SHAゲート成功 #44666182807](https://github.com/zuizui0223/alan/actions/runs/44666182807) • [完全ソース結果](../results/all_species_forest/biotic_interactions/anic2022_13species_Kew_PLANTATT_F_source_verified_20261120.json) • [解析コード](../scripts/analysis/alan_anic2022_13species_kew_plantatt_20261120.py)。

## B. Giavi, Fontaine & Knop 2021 — 昼間訪花

同じ2016草地計画ではない、独立したスイスの**6明暗草地ペア**実験（Nat Commun DOI 10.1038/s41467-021-22011-8）。元Zenodo DOI 10.5281/zenodo.4540407を含む原21植物名と種別の照明−非照明のlog(1+訪花数)回帰係数を、main保全済み元集計（Git blob SHA1 b54cadf945edf2eac129d1b45ce9d5edcbda78cf）から取得。

属止まりsp、aggr、s.l.を除き真正二名法は10植物種。独立Kew BIFloraExplorerからPLANTATT Fが8種、別ソースDöring Fが9種で厳密に照合。種間Fと訪花照明係数について、同じ種にF値ラベルを全数置換した**数学的帰無**（PLANTATT 420通り、Döring 15120通り）を比較した。**これは元6草地の照明ラベルを置換したランダム化検定ではない**。

| F形質 | 照合可能植物 | Spearman rho | 種ラベルの2側厳密置換 p |
|---|---:|---:|---:|
| 英国 PLANTATT F | **8種** | **−0.7715** | **.0381** |
| 英国 Döring F | **9種** | **−0.8806** | **.00331** |

これは「湿潤性が高い種ほど、同一照明草地実験の昼間訪花が低い」という仮説と整合的。ただし種別元係数は観測地点数が著しく異なる。

**反復地点条件による厳格な感度分析：** 照明・非照明の両方で**2以上の同一草地ペア**があった植物だけにすると、PLANTATTは4種、rho−0.7379、厳密p=.3333。Döringは5種、rho−0.7182、厳密p=.1667。向きは残るが推測精度・外挿支持は不足。

さらにZenodoの原データRDataを**著者公開バイトMD5 f60b21d82c163cbb53520fda0bc61985**で照合し、そのDATテーブルを本研究で直接読み込んだ。元は**1,884訪花行、21植物タクサ、6独立明暗草地ペア、37草地×反復測定単位、3昆虫ギルド（ハチ目、ハエ目、甲虫目）**。真正報告の訪花ゼロは**1,441行**で保持、未報告の植物・ギルドをゼロとはしていない。元測定に花量plabがあり、これを種別比較の際に調整すべき。原草地×反復×植物で両光処理が実際にある**193件**、植物×昆虫群では579件。だが**湿潤種F>7と低湿潤種F≤7が同じ草地で照明両区に共存するのは6草地中2ペア**。各草地で湿潤2種と比較2種まで揃うのはPLANTATT **1草地**、Döring **2草地**のみ。したがって群集/草地内の湿潤×照明交互作用を十分な照明独立反復で推定できない。これは本来の帰無を棄却しなかったという意味ではなく、識別不能というSOURCE-SUPPORT STOPである。

[種別係数と形質の実行結果 #44667088538](https://github.com/zuizui0223/alan/actions/runs/44667088538) • [保存した完全ソースF比較](../results/all_species_forest/biotic_interactions/giavi2021_original6meadow_kew_PLANTATT_moisture_response_20261120.json) • [真正元RData実行 #44667379311](https://github.com/zuizui0223/alan/actions/runs/44667379311) • [元RData 6草地共通支持](../results/all_species_forest/biotic_interactions/giavi2021_original_Zenodo_1884rows_6meadow_true_source_support_20261120.json) • [真正原本読取](../scripts/analysis/alan_giavi2021_originalRdata_meadowpair_support_20261120.py)。

## 論文Resultsにして良い新しい部分／してはいけない部分

**Resultとしてよい：** 二つの独立ALAN照明研究の実際の植物生産・昼間訪花と、別途定義した湿潤性の**負方向の種間順位の対応**（Anicでは3比較補正p=.0766、Giaviは緩い比較p=.038/.003、厳格現地共通支持p=.333/.167）。これは本研究で追加した透明性ある独自の**探索的二実験の定量比較**である。Figとして二つの異なるendpointを別軸にし、原独立照明ペア数と共通支持地点数を明示する。

**Resultとして不可：** 「湿潤植物がALANで直接減る」「光が湿潤植物の繁殖を抑える」「送粉・水生昆虫・捕食者が湿潤偏りを媒介する」「二研究で独立した湿潤種確証的再現」「両研究を統合したメタ解析で一般則確認」。Anicの光コードの元意味が未特定、Giaviの場所レベル比較が2しかなく、元論文はそれぞれ植物バイオマスと昼間訪花を観測しただけ。実測種子・水文変化・夜行性送粉・同一地点の植物減少は測定されていない。

### 次の機構的検定に必要な実験条件

同一原実験で4種以上の湿潤F>7、4種以上の比較F≤7を含み、**異なる照明設置または草地で5地点以上**湿潤・比較が共存するsource-originの草地実験を最優先する。単独のF種間比較で広域の種間相関を報告するだけでは元の湿潤偏りの因果機構を同定できない。花量、原観測努力、照明濃度、種・系統、実験区・期間をモデル化し、訪花数、同種花粉送達、種子・結実を分離する。これを満たさない場合は未同定と宣言し、元の主結果を質的な群集偏りに限定する。

**Disposition: DRAFT only. do not modify original alan/main manuscript or promote any causal mechanism.**
