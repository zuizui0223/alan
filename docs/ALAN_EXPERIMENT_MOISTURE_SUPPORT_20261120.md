# 実験効果と植物湿潤性の結合：現時点の支持範囲

2026-11-20。マクロ主結果で使った著者配布EIV表を同じSHA256で照合し、原論文で種が特定できる植物応答に完全二名法一致で結合した。[Tichýら原論文](https://doi.org/10.1111/jvs.13168)、[著者配布CSV](https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv)。最新版への置換や類義名の自動補完は行っていない。

| 原資料の学名 | MOIST | 実験応答 |
|---|---:|---|
| Asclepias syriaca | 5.1 | Heyの最終生物量、水追加×密度4対比 |
| Hydrilla verticillata | 12.0 | Xueの単植最終生物量 |
| Myriophyllum spicatum | 11.9 | Xueの単植最終生物量 |
| Lotus pedunculatus | 7.9 | Bennieの最終花序数。ただし暗区平均ゼロでlnRR未定義 |

Myriophyllum aquaticum、Vallisneria natans、Psychotria asiaticaはこの完全一致で値が得られなかった。Aporosaは原データと最終論文で種名が食い違うため、種特性の結合を保留した。欧州EIVの欠測を低い湿潤性やゼロと解釈しない。

有限の効果量とMOISTがともに得られるのは6行・3種・2実験プログラム。Heyの4対比は同じ1種であり4独立種ではない。Xueの2種は同じ照明ゾーンを共有し、MOISTも11.9と12.0で差が小さい。現状の高湿潤性と水生容器実験は重なり、MOIST効果を実験設定や種差から分離する支持が不足するため、湿潤性メタ回帰は実行しない。

群集総生物量や複数植物への訪花合計に、任意の代表植物のMOISTを付与していない。種の湿潤性、実験内の冠水・水追加、座標周辺の水面はそれぞれ別の説明変数である。後者の衛星由来調整変数はこの結合では未算出。

結果表は`site_context/plant_species_moisture.csv`、効果量との対応は`site_context/plant_moisture_effect_support.csv`。これは機構が否定された結果ではなく、現時点で何を独立に検定できるかを示す支持範囲の監査である。
