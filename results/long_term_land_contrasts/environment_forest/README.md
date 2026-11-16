# 長期増光・気候・市街化の種別予測比較

対象は既存の長期傾き推定が可能な427種。応答は座標ごとに同じ暦月の実調査5年以上・期間20年以上から推定した植物記録率の線形傾き。衛星夜間光の傾きも同一系列から計算する。種を選ぶ際には予測値に加え、既存の区画内線形推定の符号と、増光地点での平均記録率変化の符号を照合する。これらは同じ資料に由来し独立な再現ではない。

## 環境データ

NOAA PSL が配布する NCEP/NCAR Reanalysis 1 の地表気温と降水率。1992–2024の各地点の調査年・暦月に最近傍グリッドを対応させ、傾きと期間平均を同じ月・情報源順序で座標へ集約。気温格子は2.5度、降水は約1.9度で小島の局地気候を解像しない。降水量の単位は日換算mm、傾きは10年あたり。観測年内の気候差・非線形な履歴のすべてを入力するモデルではない。

- [気温月値](https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc) SHA256 57FFDC5F0B57FDE8ABE28A838B252455307AA0C6DE0846B8F3233FB44DE57493
- [降水率月値](https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/prate.sfc.mon.mean.nc) SHA256 0DD5BF2985DB448C47B73DF2FE0C8000FD75638A9B4C9A91F6A8D905C2517780
- European Commission JRC GHS-BUILT-S R2023A、Mollweide 1 km。衛星に基づく構築物被覆は、[1990年](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E1990_GLOBE_R2023A_54009_1000/V1-0/GHS_BUILT_S_E1990_GLOBE_R2023A_54009_1000_V1_0.zip) と [2020年](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000/V1-0/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0.zip) の差を3で割って10年変化に換算。森林伐採、農法、草刈り等を表さない。入力zip SHA256 は順に 3F2BFD83B0B84427B5521AFE98245530D0822EC7B70C1582C110326E2C41EC7E および D61CE2E88CF7441B062C0FBD9624C98618C682F4F6BEE5522663669D4CD11FE0。
- 光は既存の調和済みDMSP/VIIRSの0–63 DN。傾き・開始時の光は既存解析と共有。NDVI・開花・蛾はこのRFに含まれない。

## 空間予測

種ごとに座標を100 km区画で5分割し、未使用の区画で検証。夜間光・気温・降水・構築物被覆の4つの変化指標を含む森林と、それぞれ1指標を外して再学習した森林を同じテスト区画で比較する。全モデルには初期の光・構築物被覆、期間平均気候、調査植物種数の変化、観測初年・幅・年数、島区分も入る。評価は各テスト区画を等重みとする。skill は訓練側平均値で予測する基準に対するMSE改善割合。gain_delta は光の変化を外したときのMSE増加を同じ基準MSEで割った量。他の3指標にも同じ定義を用いる。順位はあくまで採用した4変化指標間の条件付き予測寄与。相関する変数間で寄与が分散しうる。

探索候補は、全モデルが平均値予測を上回る、光の寄与が4指標中最大かつ正、5分割中4以上で光を加えると改善、既存の光関連傾きと増光地点の平均変化がともに負、の全条件とした。追加の再分割はこの候補で同じ区画の別割付と乱数を使い96本の森林で実施。テスト地点で訓練光量範囲内の局所予測差（訓練分位点の隣接区間、他の変数を固定）も記録。再分割は独立した地域・時期の再観測ではない。予測寄与の比較は因果要因の順位・将来予測・多重検定の有意種リストではない。

再実行は4つの原資料を work/gbif_annual/ncep に置き、scripts/analysis/alan_long_term_environment_forest.py に work と results を渡し、続いて同じ引数に --repeat を追加する。
