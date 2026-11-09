# ALAN: 全球93,209,351枚からGloBI蛾訪花陽性植物の元写真ID全数回収

2026年11月9日。成果は原本写真メタデータの全件照合と抽出であり、追加で3,435万枚をダウンロードして実測花色ROIを測定したことを意味しない。

## 本当に完了した実測の出典

- [Global species-signed source count run #42607908503](https://github.com/zuizui0223/alan/actions/runs/42607908503): SUCCESS。2026-09-27版42,111名目植物分類群・93,209,351 distinct original photo IDsを、WCVP accepted botanical nameとGloBI 61,967花相互作用study recordsに一致結合。
- [16 original groups and 224 tier tar extractions run #42608966731](https://github.com/zuizui0223/alan/actions/runs/42608966731): 16/16 original source group workers SUCCESS、それぞれ14 tier separate original photo ID CSV.gz and SHA provenance receiptを保存。最後のsummary collectorだけPython import漏れで失敗したが、原本photo ID抽出データは損失していない。
- [Collector-only recovery run #42609768617](https://github.com/zuizui0223/alan/actions/runs/42609768617): SUCCESS。すでに成功した16 groupの小型receiptのみ再取得し、重い原本tarsや画像を再取得せず、**93,209,351/93,209,351**元写真IDを全数照合、**34,346,428/34,346,428**蛾正記録元写真IDを復元したことを再検証。224/224 tier source。元画像は1枚も取得せず、ROI9測色成功として0を保持。

## 既知の蛾の花訪問群への確実な全原本写真ID割り当て

| GloBI / 文献上の正記録の範囲 | WCVP accepted species | 元写真ID |
|---|---:|---:|
| GloBIに蛾の花訪問・送粉関係が記録された植物 | 1,836 | 34,346,428 |
| GloBIに蛾の pollinates 関係が記録された植物 | 493 | 13,262,462 |
| Walton 2020原著の夜間採集蛾による花粉運搬関連植物 | 41 | 3,738,985 |

ただし、GloBI蛾陽性写真のうち**32,458,089枚**はチョウ陽性記録もある植物の写真であり、蛾陽性だがチョウ正記録がない植物は**1,888,339枚**。後者も「実際に蛾専媒」とは立証されない。GloBIに記録がないことは「蛾の訪花がない」の証拠ではない。原著Walton植物も、現地の柱頭到達・種子生産や世界全域での蛾による専属受粉は未検証。

## 実際に処理するならどのくらい？

署名つき全224 tier×16 groupの元写真ID抽出結果から、GloBI蛾陽性の元写真34,346,428枚のうち**20,840,036枚はALANの同年衛星共変量に利用している原観察年2019/2021/2023/2024/2025**。全14 tier群におけるdistinct original observation UUIDの合計は**20,858,255**。同じ観察の複数写真を生物学的独立調査として重複カウントしない。

一貫した凍結済FCP YOLO11/EfficientSAM original flower ROI9方式で元写真3,435万枚すべてを直接測色するには、512枚ずつでも**67,083ワーカー群**が必要。これは実写ダウンロード・ROI算出・技術unknownの監査を伴い、メタデータの抽出とはコストが桁違い。既存のtier08 24,576 genuine original RGB acquisition run #42601821195 と3段階の回帰 run #42602113979が進行中であり、同一原本の重複取得を避ける。

## 蛾陽性植物だけでの生態系統計はすでに計算済み

[Source-authenticated old+new original ROI9 moth-restricted run #42606223917](https://github.com/zuizui0223/alan/actions/runs/42606223917) SUCCESS。旧新写真ID重複除去87,852枚の中で、蛾陽性植物の元写真試行26,420件、ROI9実測成功11,801件、厳格な花色site×year 684件、58 GSHHG島、島内異なる明暗siteを含む30島・652 site×years。島均等・原撮影年VIIRS・島/年FE・気候・人口・島面積・本土距離調整 **β=+0.66179368、nominal 95% CI [−2.35404,+3.67763]、p=.6569**。GloBI DB pollinates群は空間比較可能9島、夜間花粉運搬Walton群は0島で交互作用は推定不可。

**主結論**：写真の元IDは全数取り切れたが、現在のボトルネックは原本写真の数だけでなく、島内・同年・種組成を揃えた実測ROIと独立の夜間送粉証拠。写真に映った花領域の白色ROI9は、植物群集の白花頻度やALANによる遺伝的花色進化を直接測らない。

## 今後の取得の機械的な正しい順序

原本photo ID全数抽出 → 元のobservation UUIDとoriginal yearで1観察単位を選ぶ → 2019/2021/2023/2024/2025・GSHHG島・同年VIIRSで有効な島内異なる明暗siteへ優先配分 → taxon/size/area階級で偏りを監査 → 事前に原本source IDとsampling strategyを凍結 → 真のRGB・ROI9を分割実測 → 群別の実測島数とβ/CI・共通地域サポートを比較。効果のp値が有意になっただけで停止しない。

すでに復元された34,346,428件のsource photo IDはGitHub Actions artifactで**16群×14 gzip CSV**として利用可能。原本写真ピクセルの無条件転載ではなく、original public index内の写真ID/UUID/メタデータの派生ソース。保存期限が来ても署名つきMOIS原本source releaseと本PR抽出コードから再現できる。
