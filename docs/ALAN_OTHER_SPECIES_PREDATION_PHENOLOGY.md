# 花粉を運ぶ蛾・捕食者・開花日を植物種ごとに照合する

## 同じ一次調査から3種の蛾体表花粉を追加確認

Syskineらのコロラド州Gunnison郡4地点の2021年夏の調査で、前回の個別検証で取得できなかった蛾個体別のワークブックをDryadの固定版329747から入手した。SHA-256はDryad登録値と一致した。原表の捕獲個体番号は276件すべて別で、再捕は0件。研究者は蛾の胸、頭、伸ばした吻と小顎鬚をフクシンゼリーで拭い、採集標本・現地の花粉ライブラリーを参照して顕微鏡で形態を同定した。蛾1個体から同じ形態型が5粒以上検出されたものが記録対象。個票は植物種ごとの個体あたり粒数を含み、胸と吻のどちらに付いたかは区別できない。脚上の花粉や柱頭に届いた数は記録されていない。[原著の採集・同定方法](https://doi.org/10.1111/een.13399)・[原データとデータ辞書](https://doi.org/10.5061/dryad.j6q573npc)・[個票の公式チェックサム](https://datadryad.org/api/v2/files/3652358)

原表の花粉列名と544種の受容名が文字列で完全一致し、陽性個体がいたのは次の3種。種に同定されなかった *Artemisia* spp.、*Stellaria*／*Ribes* 判別保留の列、科名・番号だけの形態型は割り振っていない。元の早期・中期・後期の分母は91、90、95個体で、表はその期間に該当花粉を運んだ頭数。調査は灯火トラップで蛾を採集しており、灯火の有無で花に来た頻度を比較した値ではない。[季節別の個票集計](../results/shared_mechanism_crosscheck/colorado_pollen.csv)・[原表の蛾名別内訳](../results/shared_mechanism_crosscheck/colorado_pollen_carriers.csv)

| 植物 | 陽性個体（早／中／後） | 全粒数 | 元の島−本土の寄与 | 元の両側に写真がある地域年 | v0.11ラベル |
|---|---:|---:|---:|---:|---|
| *Achillea millefolium* | 0／0／1 | 11 | −0.0490 | 4 | 保持 |
| *Heracleum maximum* | 0／6／12 | 3,623 | −0.0189 | 2 | 保持 |
| *Linaria vulgaris* | 2／2／1 | 150 | +0.1543 | 0 | 旧のみ |

*Achillea* の11粒はヤガ科 *Euxoa castanea* 1頭。*Linaria* の5頭は5種で、*Sphinx vashti* 1頭の88粒を含む。*Heracleum* では中期に1,500粒、後期に2,123粒、計18頭・13種類の蛾ラベル（うち1種類は種未同定の「Noctuid spp 2」）。3種の種間で運搬の季節と担い手が異なる。蛾に依存する繁殖割合や当地での通年の主要送粉者順位に読み替えることはできない。RMBLの調査地域（38°57′N、106°59′W）は5°区画14_25で元の78地域には含まれず、2021年という暦年のみ共通する。[元の195地域年](../results/bilateral_visitors/regional_frame.csv)・[個票再集計コード](../scripts/analysis/alan_shared_colorado_moth_pollen.py)

3種はいずれも既存の544種表では運搬未確認、花粉制限実験該当件数は0。既存の確認49種と重複しないため、この追加原著を採用すれば**少なくとも52／544種**に直接証拠がある。旧台帳の版、−1.902、64種の花粉制限の推定は更新しない。*Linaria* のように旧版のみでプラスに大きく寄与する種にも、データベース版と独立した花粉運搬の証拠がある。[既存の根拠と544種の全表](../results/species_interpretation/README.md)

## 捕食者が有効な送粉を補い得る植物

*H. maximum* は共通544種に入る（WCVP 2846060、元の島5／本土17陽性地点年）。島で観察された捕食者のうち、採食の原記録を見直した後に「成虫の蛾を捕る捕食者とその蛾が同じ植物を使う」と残るのは *Dolichovespula maculata* → *Cisseps fulvicollis* → *H. maximum*。ハチ自身の記録は花上の網採集、蛾の訪花記録は2026年の日中で、上のコロラド個票には *C. fulvicollis* が0頭だった。したがって新確認の3,623粒をこの捕食された蛾が運んだとはいえない。年と陸塊を合わせて3者の所在が一致した元の地点年は1°で島1／本土8、0.05°では両側0。捕食圧や補償量を求めるには、同じ開花期間の襲撃・訪花の分母と、捕食者と蛾が届けた有効花粉の比較が必要である。[捕食の原記録](https://www.inaturalist.org/observations/32443214)・[訪花・採集の原資料](ALAN_HERACLEUM_POLLEN_COMPENSATION_20261111.md)・[空間一致](../results/dual_role_80/summary.json)

## 1年の同所と開花日の同所の違い

元の195地域年内で外部の開花・未開花調査が島・本土の明暗4区分にそろったのは共通544種中5種。これら5種各々で、同一の0.05°格子・陸塊・年の動物記録照合を再実行した。ここで対象にしたのはヤモリ2種、カマキリ、*Argiope bruennichi*、*Vespa velutina* の現地記録で、その植物上で捕食を観察した件数ではない。外部の開花調査地点は元の写真格子とずれることがある。元の格子まで限定した場合は候補4区分がそろう種が0。[再集計コード](../scripts/analysis/alan_monitoring_predator_by_species.py)・[種別照合表](../results/shared_mechanism_crosscheck/flower_predator.csv)・[厳密な格子照合](../results/regional_repeated_phenology/exact_cells/regional_counts.json)

| 種 | 状態観察件数／うち開花 | 開花と捕食者の記録が同年の同じ格子・陸塊にある地点年 | 確認された捕食者 | 開花を確認した当日の一致 |
|---|---:|---:|---|---:|
| *Holodiscus discolor* | 277／25 | 0 | 該当なし | 0 |
| *Ilex verticillata* | 99／1 | 0 | 該当なし | 0 |
| *Primula veris* | 59／53 | 5 | *A. bruennichi* | 0 |
| *Tanacetum vulgare* | 10／10 | 1 | *A. bruennichi* | 0 |
| *Tussilago farfara* | 69／54 | 4 | *A. bruennichi*、*V. velutina* | 0 |

同じ地点年が複数種に数えられるので3種の結果を足して独立な10地点年とすることはできない。全体で7地点年。春の *P. veris* や *T. farfara* では、年の一致から当日の一致に絞ると全件が外れる。前後に未開花を記録して挟めた開花期間内に捕食者が観察された例も0だった。*H. discolor* は4区分の各々で同一地点の複数日訪問がある唯一の種だが、島暗所2回はいずれも3〜4月に未開花で、他区分で開花が確認された5〜7月に同じ場所の観察を欠く。元の植物写真の島・本土双方で同地域の同地点に複数年写った3種の一つでもあるが、両側で比較できる同一の2年はない。これらの違いから、季節のずれと捕食の機会を結ぶ際には、年内の場所の一致だけで結論を出せない。[開花支持・観察努力](../results/regional_repeated_phenology/README.md)・[地域年合算の検算](../results/regional_repeated_phenology/predator_support.json)・[年別写真・衛星の結果](ALAN_TEMPORAL_EVIDENCE_SYNTHESIS.md)

再実行は `python scripts/analysis/alan_shared_colorado_moth_pollen.py CACHE.zip results/shared_mechanism_crosscheck` と `python scripts/analysis/alan_monitoring_predator_by_species.py WORK_DIRECTORY results/shared_mechanism_crosscheck/flower_predator.csv`。前者は公開ワークブックを検証し、後者は既存の反復観察と日付確認済みの動物データを利用する。
