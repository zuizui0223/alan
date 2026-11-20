# 抽出済み7研究の一次資料による座標・環境監査

2026-11-20。`effect_catalogue.csv` の7 study_id / programme_idを保持し、`database/alan_experiments/site_context/primary_coordinates.csv` に21地点行を作成した。数値座標あり16行、欠測5行。町中心のジオコーディング、地図からの座標推定、衛星抽出は行っていない。

## 公開済みの地点別座標を確認できた研究

**Giavi2021：12草地すべて。** [原著](https://doi.org/10.1038/s41467-021-22011-8)のSupplementary Information Table1（MOESM1、p2）は6組の明暗草地を緯経度で列挙する。取得した原著Source Data `41467_2021_22011_MOESM6_ESM.xlsx` の `Table S1` にも同じ12行があり、CSVにはこちらの数値を転記した。数値は4–5桁の小数度で、位置精度そのものは記載されない。研究単位は100m調査線を含む草地であり、個々の花や線端点の座標ではない。

原ファイルは `../gbif_annual/alan_experiment_sources/giavi2021/` に保存。Source Data SHA256は `9b5a77555433dac6d20f353a5cbca545ef3942f2cf6516302835c6a1807148ca`。取得は[EuropePMC公式補遺配布](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7966740/supplementaryFiles)。MOESM5 reporting summaryで採集期間2016-06-01〜09-06を確認した。本文の「未管理草地」は水位や湿潤性の測定ではない。Knopとのprogramme重複フラグは解除しない。

**香港は対照的に欠測のまま。** [Guo2025原著](https://doi.org/10.1007/s00114-025-02031-w)のFig2aには7組の位置が星印で示されるが、数値座標の表は確認できなかった。[Figshare原データv1](https://doi.org/10.6084/m9.figshare.27273918.v1)の全CSV見出し、README、Rファイルと履歴を確認したが緯経度列はない。Lamma Island、Kadoorie Centre、Little Hawaii Falls、Lam Tsuenの4名称で行を保持し、lat/lonは空欄。地名の滝・島から水辺のプロットと決めつけない。各100m²の14プロットのGPSが必要。

## 代表座標を持つ4研究

| 研究 | 原資料の座標 | 位置づけ |
|---|---|---|
| Parkinson2020 | 42°38′45.32″N, 83°33′13.48″ | Haven Hill Lakeの研究地点代表。原文で経度のWが欠落しており、明記されたMichigan USAに従い負号を付したことを記録。20メソコスム個別座標ではない。 |
| Anic2022 | 50°10′N, 5°7′W | Penryn campus代表点。砂充填の排水付き人工メソコスムであり、周辺の地表湿潤性が処理槽内水分を表すとは限らない。 |
| Lu2026 | 33°40′N, 101°52′E | Gannan研究所代表点。5組の25m²区画の位置は不明。水源涵養地域という広域名称を湿地分類に用いない。 |
| xue2026 | 29°30′N, 112°48′E | Dongting実験場代表点。6照明ゾーンのGPSではない。原著の別座標29°49′N,113°27′Eは材料採集地Jingzhouなので採用しない。 |

出典：[Parkinson Methods Study site](https://doi.org/10.1371/journal.pone.0240138)、[Anic Methods2.1](https://doi.org/10.1002/ece3.9157)、[Lu Methods2.1–2.2](https://doi.org/10.1002/ece3.74368)、[Xue Methods Study system/Experiment1](https://doi.org/10.1002/ecy.70458)。度分秒は数値変換しただけで、追加の空間精度は付与していない。分単位座標の四捨五入誤差は緯度方向だけでも最大約0.9kmになり得るため、数十mの水面近接性を識別できない。

## 座標が確認できないHey2020

[原著](https://doi.org/10.1007/s00442-020-04679-9)、補足図S2、原Excel metadataはBlandy Experimental Farmと10区画の設計を示すが、監査した資料には数値座標がない。所属住所・農場中心点を実験座標に置換せず空欄とした。2017年4週間の処理、水追加3.8L/週対自然降雨、商用培養土を記録した。自然降雨群を乾燥地と呼ばない。7月が多雨だったことから実験の正確な開始終了日を作らない。

## 共変量として使う際の境界

植物EIV MOISTは**種の環境選好**、座標からの水面距離は**実験地点の景観条件**、水槽・冠水処理は**実験内の水分条件**であり、相互に代入しない。今回の資料が直接確認する淡水環境はParkinsonの浅い湖岸（水深0.30–0.45m）とXueの水深28cmの容器である。その他の草地、季節性の雨、森林について「湿潤系」と自動判定していない。

次段階で細かい水面距離を計算するなら、まずGiaviの12地点を対象に、100m調査範囲を含む位置不確実性を保持する。分単位代表点の研究は広域気候などとの結合候補であり、正確な灯火周辺の水面距離には追加GPSが必要。香港・Heyの欠測を埋めるには著者提供の座標が必要。人工槽の水分条件は衛星の水面検出でなく一次Methodsから保持する。

7研究・21行・16数値座標を検証した。研究や独立反復が21あるという意味ではない。主結果・効果量・本体マスターは変更していない。
