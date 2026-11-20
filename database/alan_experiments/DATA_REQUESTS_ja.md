# 統合を進めるための原データ照会案（2026-11-21）

状態：未送信。公開資料を先に調べたうえでの照会文案。著者への連絡は未承認であり、このファイル作成によって送信したことにはならない。連絡先・署名は送信時に確認する。

目的は有意な結果を増やすことではなく、正・負・不明瞭な結果を同じ採択条件で評価できるようにすること。以下の回答が得られても、湿潤／乾燥の反復支持範囲が自動的に成立するわけではない。データを得られない研究をゼロ効果や不適格として扱わない。

|優先資料|統合を止めている具体的な不足|最低限必要な回答|回答後に可能になること|
|---|---|---|---|
|Lu ECE2026、10.1002/ece3.74368|訪花総数と5群合計が昼20行・夜3行で不一致。夜のセッションキー重複。JPE別報との重複不明|修正版または各列・行の正しい解釈、独立区画とセッションID、対象植物・観察努力、別報の年／区画対応|昼夜の比較を再構成し、Giaviとの対象範囲の違いを明記した採否を判断。種子形成とは統合しない|
|Macgregor2019、10.1002/ecs2.2550|正方向の種子形成成功結果について数値係数・分散を未回収。植物／区画の再利用構造が必要|解析用データとID、または条件が明示された係数・SE・対比共分散、Appendix S1|同じ繁殖エンドポイントの実験と比較可能か判断。照明主効果と照明×訪花許可交互作用を混ぜない|
|Parkinson2026、10.1111/fwb.70236|既存2023単夜データとの重複と固定トラップの7夜反復|トラップ×日付×処理×生息場所×分類群の計数・体サイズ／重量推定、June3の対応|既存系列の時間拡張。新しい独立研究として追加しない|

## 1. Lu：訪花記録と関連論文の対応

Subject: Clarification of published visitation data and experimental overlap (10.1002/ece3.74368)

Dear Dr Lu and colleagues,

I am compiling an independently extracted synthesis of field lighting experiments, separating visitation, pollen delivery and plant reproduction. Thank you for providing the workbook and analysis code for your Ecology and Evolution article (10.1002/ece3.74368).

In the supplied Date.xlsx visitation sheet, the abundance column differs from the sum of the five insect-group columns in 20 daytime rows and three nighttime rows. Three nighttime session keys (sampling_time, treatment, block, time, repeat) are repeated with different observations. I have retained these records without correcting or deleting them. Could you clarify the intended interpretation, or share a corrected version with stable plot and observation IDs? A row-level discrepancy table is available if helpful.

Could you also confirm the plants and observation effort represented by each session, and whether the plots or observations overlap with your Journal of Plant Ecology article (10.1093/jpe/rtag160)? For the latter, an accepted manuscript and the data underlying pollen and seed-set analyses would help clarify treatment replication and model scales; I have located the supplementary tables.

I will preserve all result directions and avoid counting related reports as independent experiments. Please let me know any citation or reuse conditions.

Thank you for considering this request.

## 2. Macgregor：正方向の繁殖応答を同じ基準で評価

Subject: Reproductive data and model contrasts for 10.1002/ecs2.2550

Dear Dr Macgregor and colleagues,

I am compiling an independently extracted synthesis of artificial-light field experiments and would like to represent your 2019 Ecosphere study (10.1002/ecs2.2550), including its reported positive full-night response in seed-setting success.

I located the article through White Rose Research Online and the author analysis scripts at Zenodo (10.5281/zenodo.1743917), but have not recovered the underlying reproductive dataset or Appendix S1. Could you share these, including lighting and pollinator-access treatments, observation round, plot and reused-plant identifiers, seed-setting outcomes, and seed counts where recorded?

If individual records cannot be shared, a lighting contrast with its scale, treatment reference, SE and covariance with other reported contrasts would be useful, together with the independent allocation and repeated-observation structure. I will keep seed-setting success separate from seed counts conditional on successful undamaged capsules and will not pool lighting main effects with visitor-access interactions.

Please let me know any citation or reuse conditions. Thank you for considering this request.

## 3. Parkinson：同一系列の季節的な拡張

Subject: Seasonal trap data and overlap with the June 2022 study (10.1111/fwb.70236)

Dear Dr Parkinson and colleagues,

I am compiling an independently extracted synthesis of field lighting experiments and am interested in your seasonal aquatic–terrestrial insect study (10.1111/fwb.70236).

Could you share the trap-level observations across dates, with stable station IDs, light treatment, habitat placement, sampling date and duration, taxon counts, and the measurements or conversion factors used for biomass? A data dictionary and missing-sample flags would also be valuable.

The site, treatment design and June 3 sampling appear to overlap with Parkinson and Tiegs (2023; 10.1098/rstb.2022.0361), whose public workbook I have examined. Could you confirm the overlap and how station IDs correspond? I would treat repeated dates as observations within an experimental programme, rather than additional independent studies.

I will distinguish trap catches from insect production and preserve all response directions. Please let me know any citation or reuse conditions. Thank you for considering this request.

## 回答後の処理

受領資料は公開資料と分け、利用条件・バージョン・受領日を記録する。旧データを上書きせず、修正行・ID対応・原著結果との一致を検証する。独立照明区の平均と共分散を復元できた段階で効果量を抽出する。個別記録が非公開の場合は、許諾された派生集計・再現コードのみ公開する。

モデル対比が提供された場合は、尺度・対比の定義・独立単位・共分散を確認し、群平均lnRRとは別に保持する。比較可能性を確認できれば同じ尺度の対比として統合を検討できるが、SEのみを使って対応構造を無視しない。

現時点では照会送信数0、回答0、新規メタ効果0。これは依頼内容を具体化した文書であり、データ回収・採択・メタ解析の完了を表さない。
