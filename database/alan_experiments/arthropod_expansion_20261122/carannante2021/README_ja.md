# Carannante et al.（2021）：4独立河川明暗トラップペア原データの取得ゲート

2026-11-24の**新しい独立実験候補の一次本文確認および公開原データ取得監査**。原著：[Carannante, Blumenstein, Hale & Arlettaz（2021）, *Ecological Solutions and Evidence*, DOI 10.1002/2688-8319.12053](https://doi.org/10.1002/2688-8319.12053)。研究名は *LED lighting threatens adult aquatic insects: Impact magnitude and distance thresholds*。作者本文のData availabilityが[University of Bern BORIS record 151803](https://boris.unibe.ch/151803/)を指定する。

## 物理実験設計を原著から特定

- スイス・Gürbe川の2017年夏。**Study A**は4つの照明＋暗トラップの対応ペア。ペア同士は200 m以上離れ、同一ペア内の明暗は90 m離れていた。2ペアは水際から3 m、2ペアは40 m。これを**6夜反復**した。真の明暗ペア数は4、**48のトラップ×夜イベント、昆虫約10万個体は独立した照明実験反復ではない**。
- 20 W／1500 lumen白LEDを夜間約2.5時間使用。水生幼虫をもつと考えられるカゲロウ目、トビケラ目、体長3～8 mmのハエ目を対象に、飛翔インターセプト＋粘着トラップ（前二目は手取りも）を組み合わせた。**Dipteraを完全に水生由来と同定していない**。真の羽化トラップでも、植物への訪花者の確認でもない。
- Study Aの原著記載捕獲数は合計99,852個体（Ephemeroptera **181**、Trichoptera **2,146**、Diptera **97,525**＝**約97.67%**）。著者の概略一時間当たり平均は明区412個体・暗区4個体。しかし種群・努力・採集夜の分散を回収していないため**412/4を新しい正式メタRR／SEへ変換しない**。
- **Study B**は照明トラップの河川距離3、10、20、40、60、80 mでの勾配試験であり、**暗区対照を持たない**。A/Bで同じ河川・研究グループを共有する可能性があり、単純に2独立実験とはしない。距離勾配は「照明／暗区の距離×照明交互作用」と同義ではない。

## 今回取得できた出典と新しいアクセス境界

公開BORIS旧URL `https://boris.unibe.ch/151803/` をGitHub Actionsの公共アクセス環境から開くとHTTP200で**https://boris-portal.unibe.ch/items/f559ab1e-e693-4fbc-b7fe-683134f6dbe1** へ転送された。**University of Bern BORISの実際の移行先UUIDが判明**した点は新しい取得進展。

ただし、正式の公開REST検索2ルートと大学ポータル検索は**HTTP200でも同一のHTMLレスポンス（各7,308バイト）**になり、機械可読JSON・採集値・元のCSVファイルを取得できなかった。BORISは匿名アクセスを制限する場合があり、HTTP200を「データ取得成功」と誤認してはいけない。今回の専用ゲートはAPIにHTMLが返った場合を`API_RETURNED_HTML_SHELL_NOT_AUTHOR_DATA`、アクセス制限の文言を検出した場合を`PUBLIC_ACCESS_RESTRICTED_HTML_SOFT_200`として区別する。**アクセス制限を迂回せず、実際の4ペア×6夜原データ未取得**。脚注とMethods/Supplementの項目があることだけで、元データが計算可能だとすることはできない。

[公表済み記事](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1002/2688-8319.12053)、[大学BORISの公開移行先](https://boris-portal.unibe.ch/items/f559ab1e-e693-4fbc-b7fe-683134f6dbe1)。

## メタ解析採択までの客観的ゲート

1. 大学の正規公開経路または許可を得た著者配布で**Study Aの4ペア×明暗×6夜×3目**の原データ・採集努力・照明距離・原ファイルSHA256を回収。**本段階の原データ確認はSTOP**。
2. 真の4ペアを単位に夜反復を集約し、各分類群の照明／暗比または生計数モデルとその分散を算出。独立反復夜を標本nにしない。1照明設置が6夜で独立6基に増えるわけではない。
3. カゲロウ目／トビケラ目と3–8 mmハエ目の**相対捕獲応答**を同一ペアで比較できる。ただし`aquatic origin versus terrestrial origin`という既存のParkinson捕獲メタの直接kへは使えない（ハエ目の幼生生息地未確定）。
4. 水際3 m／40 mの2ペアずつだけで照度・距離修飾を検定できても不精確。Study Bの距離勾配を明暗交互作用へ転用しない。

**状態：原著実験設計の確認済み・BORIS出典URLとUUID確認済み・原データ取得保留。新規メタ採択k=0。** [実行可能な公開アクセスゲート](../../../../scripts/acquisition/alan_carannante2021_boris_source_gate_20261124.py)／[GitHub Actions source run](https://github.com/zuizui0223/alan/actions/runs/45332689522)。既存108行と4探索的統合は不変。
