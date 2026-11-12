# 北米北東部のトウワタとアザミを実際に訪れた昆虫

対象は地域21_26（西経75～70°・北緯40～45°）と元の2019・2021・2023・2024・2025年。地域内で訪花が実際に記録され、どの陸塊・元の撮影地点年まで結び付くかを確認した。主要結果のGloBI版・分子・モデルは変更していない。

## 地域と元の地点年での直接訪花記録

|植物|5°区画での訪花記録 島／本土|うち蛾類 島／本土|昆虫の記録種数 島／本土（共通）|同じ種の元の0.05°×陸塊×年に一致する訪花記録 島／本土|
|---|---:|---:|---:|---:|
|*Asclepias syriaca*|165／1,627|19／195|50／215（39）|0／8|
|*Cirsium vulgare*|12／97|0／2|9／28（6）|0／0|

トウワタでは同じ蛾 *Eichlinia cucurbitae* が島8・本土46件、マルハナバチ *Bombus griseocollis* が39・234件記録された。*Hemaris diffinis* も両側で記録されている。アザミでは *Bombus impatiens* が2・18件で、両側の訪花者に重なりがある。地域内で両側に訪花記録がある年はトウワタ5年すべて、アザミは2025年を除く4年。元の地点年で一致したトウワタ本土8件の内訳は膜翅目5、チョウ2、蛾1（*Ctenucha virginica*）である。植物撮影日まで同じものは0件。これらは観察者による陽性記録で、撮影努力・調査時間・花数が不明。蛾類には昼行性の種も含まれる。局地的な主送粉者、夜間訪花頻度、花粉を届けた量、島での送粉者交代はこの件数から推定できない。[トウワタを訪れた蛾の島記録](https://www.inaturalist.org/observations/172327972)・[本土記録](https://www.inaturalist.org/observations/32449602)・[アザミを訪れたマルハナバチの島記録](https://www.inaturalist.org/observations/233154743)・[本土記録](https://www.inaturalist.org/observations/89831597)・[元の地点年の蛾](https://www.inaturalist.org/observations/28808103)

2026年11月12日の現行GloBI APIで各植物の地域内の visitsFlowersOf 観察行を取り出し、対象年のiNaturalist原観察ID 2,856件を再取得した。原票の「Interaction->Visited flower of」（field 3126）に植物の正確な学名があり、元の解析年、公開・非秘匿座標、位置精度0 m超1 km以下、research grade、昆虫が種以下に同定される条件を適用。島・本土の判定は元のGSHHGで行い、1,901件の植物×原観察が残った。訪花者はレコードID×植物で重複除去し、亜種名を種に集約した。新旧GloBIの陽性ラベルを再比較するための検索ではなく、日時と場所を持つ原票を探す索引に使っている。今回の切り出しは独立の調査努力をもつ標本でも、地域の全観察の網羅的台帳でもない。[GloBIの観察データAPI](https://globalbioticinteractions.org/data/)・[原票の訪花植物欄](https://www.inaturalist.org/observation_fields/3126)・[種名別・陸塊別の集計](../results/northeast_local_visitors/visitor_taxon_counts.csv)・[年・陸塊・空間条件ごとの結果](../results/northeast_local_visitors/verified_visitors.csv)・[データ照合コード](../scripts/analysis/alan_northeast_focal_visitors.py)

## 区画内の時間付き調査で確認できること

米国NRCS・ロードアイランド大学・コネティカット農業試験場の保全草地調査では、表1の26地点がすべて対象5°区画内に入る。報告書の小数2桁の座標を元のGSHHGに照合すると島3草地（GO1、AG1、HH1）と本土23草地。対象植物が撮影されている0.05°地点を各調査年に照合した一致は0で、小数2桁の丸め誤差±0.005°まで許容しても該当しなかった。報告書はトウワタを播種種、アザミを自生侵入種として扱っており、保全草地という限定された環境での比較である。[調査方法・座標・植物別結果（表1・3・5）](https://www.nrcs.usda.gov/sites/default/files/2026-02/Assessing%20Plant%20and%20Pollinator%20Response%20to%20Practices%20Installed%20Through%20the%20New%20England%20Pollinator%20Partnership%20in%20Rhode%20Island%20and%20Connecticut%2C%202025.pdf)

2023–24年の植物別10分調査（ミツバチ・マルハナバチ・*Xylocopa virginica* を除くハチ）では、トウワタで11地点15調査・150分に304個体22種を捕獲し、アザミで3地点3調査・30分に16個体7種を捕獲した。2023–25年の別プロトコルでは、トウワタへの訪花をミツバチ49、マルハナバチ132、*X. virginica* 5、アザミへの訪花をそれぞれ0、6、0個体記録した。大きな3分類群の調査時間は草地全体に対する時間であり、植物ごとの分母にできない。植物別の表は各年・全草地をプールしており、島3草地のどこでどの植物に観察時間を使ったか、訪れた昆虫がどちらの陸塊で何個体だったかを復元できない。原報には各採集容器に日付・地点・植物・開始時刻を記録し、標本を大学・試験場に寄託したとある。原報の付表、研究室のプロジェクトページ、DataCiteのプロジェクト題名と著者名による検索から、この地点別原票の公開データは確認できなかった。夜の蛾、単回訪花後の花粉付着、種子生産も対象外である。[研究室の公開ページ](https://web.uri.edu/beelab/research/assessing-pollinator-response-to-nepp/)・[研究者の土地被覆研究抄録](https://www.eaglehill.us/NENHC_2026/program/NENHC2026-poster-abstracts.pdf)

以前確認したBudburst／UGAの候補表でも、この区画・対象年には2種の局所的な時間付き記録がなかった。区画内の他植物に対する有効なBudburst調査は本土に60件ある。ニューヨーク市のGreat Pollinator Projectは2007年開始の4年の訪花調査を2023年に論文化したもので、刊行年を元の5観測年には数えない。以上から、島と本土で実際に訪れる昆虫名の共通性は示せるが、夜間光への反対方向の構成比寄与を、有効送粉者の交代で説明する局地的な比較は成立しない。[既存の調査原表](../results/bilateral_visitors/us_candidate_events.csv)・[地域別調査件数](../results/bilateral_visitors/regional_frame.csv)・[NYC訪花研究](https://doi.org/10.1016/j.landurbplan.2023.104689)・[実施機関による調査開始記録](https://www.amnh.org/research/center-for-biodiversity-conservation/research-and-conservation/biodiversity-exploration-and-monitoring/past-projects/metropolitan-area)
