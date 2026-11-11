# 蛾の捕食者と島の夜間光　固定した植物群集解析への接続

## 元の5年比較とヤモリの地理的重なり

蛾の訪花記録を持つ植物種の割合について、島−本土の夜間光に対する傾きの差 −1.902 pp / log(1＋同年VIIRS) を主結果とした。旧GloBI訪花記録、観察植物のWCVP受容種分母、5年、島933・本土9,765地点年を再構成し、全地点キー集合のSHA256と元係数を照合した。GloBI v0.11に替えると同じ地点で +0.750［−1.076, +2.577］だったという既報も解釈に含めた。[固定済み主解析](ALAN_GLOBI_PRIMARY_RESULT_FREEZE_20261111.md)・[記録更新の結果](ALAN_GLOBI_V011_SAME_FIVEYEAR_20261111.md)

PR #40の2019・2024年の同所性監査を、元の5年すべての固定地点に接続した。対象は元監査が指定した *Hemidactylus frenatus* と *H. mabouia*。GBIFの人による観察30,921件を取得し、種同定・検出状態・座標誤差が0より大きく1 km以下を満たした18,799件のうち17,852件に元解析と同じGSHHGの陸塊を割り当てられた。1度格子を共有する隣島は同じ陸塊として扱わない。同年・同一陸塊・同じ1度格子でヤモリの記録があった植物地点は島150地点年（30島）、本土1,676地点年。植物の元の0.05度格子まで一致したのは島60・本土443地点年だった。この一致にも同日・同じ花での観察という条件は含まれない。[元の2年監査](../results/alan_globi_full_v011_gecko_same_island_photo_source_20261111.json)・[GBIF occurrence API](https://techdocs.gbif.org/en/openapi/v1/occurrence)

1度格子でヤモリ記録がある地点のうち、同じ年・5度地域で島と本土の両方に比較対象があったのは島88・本土1,236地点年（34年地域、15地理ブロック）。この範囲で元解析の共変量、年地域固定効果、各年地域の島本土を均衡化する重み、空間クラスタ区間を適用した。島−本土の傾きの差は **−2.803 pp［名目95%区間 −8.803, +3.196］**。島側の傾きは −2.098、本土側は +0.705。34固定効果を明示した別のWLS実装でも係数と標準誤差（2.797）が一致した。観測努力と分布により選択された部分集合での関連であり、未記録地点をヤモリ不在とした比較は行っていない。−1.902からの点推定の差をヤモリの寄与率や捕食による増幅と解釈する根拠はない。[係数・年別地理支持・GBIF出典](../results/alan_gecko_primary_fiveyear_20261111.json)

## 蛾の捕食と花利用を分ける

GloBI v0.11全24,453,954行の厳密な集計では、蛾科の獲物への摂食・捕食関係が41,270行あった。データベースの生活段階欄に基づく分類では成虫4,420行、成虫以外547行、不明36,303行である。関係行は独立な捕食イベント数と一致するとは限らない。蛾の捕食が記録された動物種名は鳥類349、クモ類155、コウモリ94、ヤモリ7、ムカデ類3。同じ種に花訪問または送粉のDB記載があるのはそれぞれ84、30、6、1、0だった。*H. frenatus* には蛾の捕食と送粉関係のDB記載がある。*H. mabouia* には蛾の捕食が記録されている。生息地、行動時刻、獲物の発育段階、花での花粉移送の確認は別の証拠を要する。[12群・3出典区分の種別表再照合](../results/alan_strict_predator_roles_verified_20261111.json)・[原走査](https://github.com/zuizui0223/alan/actions/runs/42952223131)・[後段の項目名修正後の認証](https://github.com/zuizui0223/alan/actions/runs/42954755175)

オーストラリア北東部の野外では *H. frenatus* が在来の *Gehyra dubia* より幅広い光環境を利用する。これは照明の近くで昆虫と出会う機会が変わり得ることを示す。本研究で問題となる機構は、照明に集まる蛾のうち実際に対象植物へ花粉を運ぶ成虫への捕食が増え、その結果として植物の結実と次世代の加入が変わる、というものだ。幼虫の捕食なら植物の被食を減らす経路もあり、捕食者自身が送粉に寄与することもある。したがって、捕食者の共存だけで植物への影響の符号は決まらない。[Zozaya et al. 2015](https://doi.org/10.1111/aec.12287)

## 捕食された蛾の種名と、訪花植物の種名をつなぐ

GloBIの元ファイルをSHA256で確認し、全24,453,954行からヤモリ2種が蛾を食べた関係（eats/preysOn とその逆向き）を抽出した。*H. frenatus* は7行・種名のある獲物5種、*H. mabouia* は12行・同4種。生活段階が成虫と分類された名前付き獲物は各3種で、それぞれ1種が、−1.902に使った元の訪花記録に学名で一致した。両ヤモリの捕食19行はすべてiNaturalistまたはtrophiCHに該当する出典群に属する。独立な現地調査が19件あるという意味ではない。[実データ走査](https://github.com/zuizui0223/alan/actions/runs/42953901880)・[獲物と発育段階の原表](../results/alan_two_house_gecko_species_prey_20261111.csv)

| 捕食者と名前の一致した成虫の蛾 | 旧記録で結ばれる植物名 | 主解析で観察された受容種 | 3者が同年・同陸塊・同1度格子に記録された地点年（島／本土） | 同0.05度格子（島／本土） |
|---|---:|---:|---:|---:|
| *Hemidactylus frenatus* と *Spoladea recurvalis* | 106 | 88 | 9／69 | **2／5** |
| *Hemidactylus mabouia* と *Feltia subterranea* | 3 | 3 | 0／1 | **0／0** |

3者照合では**植物と、その植物に訪花記録のある蛾の成虫と、その蛾を食べた記録のあるヤモリ**を種名単位で限定した。蛾2種のGBIF人観察21,908件のうち精度・検出条件に合う13,889件、陸塊に割当可能な13,453件を使い、表ではGBIFにAdultまたはImagoと明記された蛾に限定している。旧訪花表で *S. recurvalis* は109行（訪花106、送粉関係ラベル3）、*F. subterranea* は訪花3行。粗い格子から植物観察の格子に絞ると、*S. recurvalis* の3者がそろうのは5植物種・7地点年である。場所と年がそろっても出来事の同時発生、同一個体の花粉輸送、照明下での捕食は未確認で、GloBIとGBIFには観察出典の重複もあり得る。この支持範囲から捕食媒介係数を推定することはできない。

## 再実行とデータの来歴

固定入力は元の5年解析のINPUTS・SIGNEDとGSHHGの海岸線を用いる。GBIFはライブの観察APIで取得したため、取得時の圧縮スナップショットのSHA256、採用したデータセットのGBIFキー、検索種・年・地理範囲を結果JSONに保存した。再取得時は件数が変わり得る。個別の観察座標は非公開の作業ディレクトリに保存する。

```sh
python scripts/acquisition/alan_gecko_gbif_fiveyear.py PRIVATE/gecko
python scripts/analysis/alan_gecko_primary_fiveyear.py INPUTS SIGNED PRIVATE/gecko gecko.json --write-blocks PRIVATE/blocks.json
python scripts/acquisition/alan_gecko_gbif_fiveyear.py PRIVATE/moths --linked-moths
python scripts/acquisition/alan_gecko_gbif_fiveyear.py PRIVATE/multiguild --multiguild --blocks PRIVATE/blocks.json
gh run download 42953901880 -R zuizui0223/alan -n two-house-gecko-moth-prey -D PREY
python scripts/analysis/alan_gecko_primary_fiveyear.py INPUTS SIGNED PRIVATE/gecko final.json --prey PREY --moths PRIVATE/moths --multiguild PRIVATE/multiguild
```
