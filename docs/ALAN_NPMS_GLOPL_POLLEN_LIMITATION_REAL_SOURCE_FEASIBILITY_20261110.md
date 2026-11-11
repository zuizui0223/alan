# ALAN: NPMS反復植生調査 × GloPL花粉制限実験の真正原本接続（2026-11-10）

## 結論

**実測花粉制限はALAN研究へ導入可能。ただし (1) NPMS調査区での植物被度・出現変化と (2) GloPL実験地点での人工授粉による繁殖改善を別の応答として扱う。** 二つのデータベースを「同一場所の花粉制限の媒介分析」とは呼ばない。

### データ粒度と原本

- NPMS 2015–2025 EIDC DOI https://doi.org/10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8。原データZIPの外側はサーバーの再生成によりSHAが変化するため、7つの内部CSVについて固定したサイズ・CRC・行数で照合。[source workflow #42753594112](https://github.com/zuizui0223/alan/actions/runs/42753594112) SUCCESS。
- NPMS植物観測行 **271,701**、サンプル情報 **27,611**、NPMS locations元表 **5,360**（5,361番目は未識別`NA`地点IDの疑似群）、属性 **293,561**。**8サンプル**が真のlocation IDを持たず、NPMS地点年の解析から除外。再利用可能な測定地点は**5,360**、有効な訪問サンプルは **27,603**、少なくとも2年で再訪した地点は **3,358**、少なくとも3年では **2,343**。
- GloPL 公開データ https://github.com/idiv-biodiversity/pollen-limitation-data-descriptor/blob/abbe193770981b2c61713e75a0f01e004bdded8e/Output/GloPL.csv（Dryad原本と同一、SHA256 `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`）。**2,969**花粉添加実験、**1,265**原本受理種名、**1,248**実験地点。PLは開放状態に対する人工花粉追加の生殖生産増加の公表log反応比であり、蛾の送粉への依存度や植物個体数そのものではない。
- GloPL植物学名はもともと`Solidago_canadensis`のようなアンダースコア表記。**アンダースコア→スペースの書式のみ正規化**、学名の異名を自動類推せず、二名法が明瞭な **GloPL 1,223種 vs NPMS 1,236種**をexact名前照合。**共通79種**。NPMS **3,920**有効調査地点、**21,032**種×地点×年記録。79種のうち、独立したGloPL研究が**2件以上ある種は19種**、このうちNPMSで20地点以上に出現したのは**11種**。GloPL未掲載種のPLを0と扱わない。

### 直接のGloPL実測PL×実験前年までのPCNL

- GloPL原本には論文の掲載年`Year`とは**別に、実測年`year_of_study`**が存在。観察バイアス／光量結果を見ずに1993〜2015年の**1,419**元実験レコード、**597地点、853地点×実験年、341研究、507植物種**を先に凍結。1992年から光量があるため、1992年実験は事前のALAN曝露を設定できず除外。より広い1992〜2015年の1,444件とは混同しない。source cohortのSHA256 `ac8e8d191edbe175fa0f72cb61387eece986d31487ae2db9cdf3619c58d80930`。[原本検証 #42753594112 artifact](https://github.com/zuizui0223/alan/actions/runs/42753594112/artifacts/13734218883)。
- 現在、Zenodo PCNL v3の実物1992〜2014年23ラスタを**公式MD5照合**し、実験年より前の光だけを6つの互いに重ならない年度群として原本座標へサンプリングする[Actions workflow](https://github.com/zuizui0223/alan/actions/runs/42754524785)を実装。出版年、実験時の地上照度、植物繁殖指標を光量に代用しない。
- 実際のPCNLが全て揃うか、研究・地理的クラスターごとに明暗勾配が存在するかのsource支援を検証したあと、**研究均等重み、論文クラスタ誤差、実験年・緯度経度・実験方法補正**で探索的なPL係数を計算する仕様を事前固定。**地理ブロック感度も併記。光量が同じか情報不足なら「推定不能」として扱う。**
- この関連を示せても、ALANがPLを**引き起こした**ことや**蛾媒植物の減少**とは言えない。実験研究の地点選択、土地利用、気候、花粉提供者集団の違いは大きな交絡。

### NPMSの植生変化へのPLモデレーションで注意

- NPMSは`Wildflower/Indicator/Inventory`等の異なる調査種の対象選定・被度記録があり、**調査水準を揃えずに未記録を種の不在としない**。同一plot ID・再訪年・調査法・春秋の季節を合わせたcover/occupancy frameが必要。
- GloPLのPL値は別時代・別地点の**種レベルの外部実験結果**。NPMSの調査区で測ったPLとは異なる。種ごとの研究間PLばらつきとNPMSの被度測定系を整備してから、ALAN×種のPL傾向のモデレーションを試す。
- NPMSのチャネル諸島は**14実際のNPMS地点、そのうち9地点が2年以上**。英国は地理的にも大きな島であり、他の本土大陸の調査区は含まれない。**全球島／本土三重交互作用はNPMS単独では識別できない**。

## ALANでの結果管理

現在のPR [#31](https://github.com/zuizui0223/alan/pull/31)は、取得と源データ品質・標本数・交絡の監査です。GloPLのPL×長期ALAN係数やNPMS植物個体群への因果影響が原本実測と合わない限り、肯定的な進化・島嶼効果を付け加えない。

使用コード：
- [NPMSとGloPLの原本形式](../scripts/acquisition/alan_npms_glopl_source_overlap_preflight_20261110.py)
- [NPMS反復調査×GloPL種のexact overlap](../scripts/analysis/alan_npms_glopl_actual_species_repeat_support_20261110.py)
- [実際のGloPL実験年による事前標本凍結](../scripts/analysis/alan_glopl_true_experiment_year_ALAN_preflight_20261110.py)
- [GloPL×PCNL年度別実原本結合](../scripts/analysis/alan_glopl_prior_PCNL_6year_shard_20261110.py)
- [GloPLの実測PL相関モデル](../scripts/analysis/alan_glopl_real_PL_priorPCNL_exploratory_20261110.py)
