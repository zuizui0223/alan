# ALAN湿潤植物の機構的Resultsを原実験13種に求められるか：Anic 2022の真正原数値再解析

2026-11-20 / Draft PR #88 / **MOIST×照明感受性についての元資料共通支持不足を実測で判定**

## 科学的な検定の位置付け

元 \`alan/main\` の主結果は、1992–2024年の長期記録3,901評価種から抽出された「増光に関連する記録率低下候補」が比較植物より湿潤生育地を好むという種群偏り（環境調整後のMOIST差 +0.534、主要地域ジャックナイフCI [+0.158,+0.910]）。これをALANの直接因果効果や花粉・種子低下に置き換えてはいけない。

このノートの独立した機構的問いは、**同じ照明操作を受けた多数の植物種について、種の湿潤性が植物バイオマスの照明応答の符号・大きさを予測するか？** である。古い群間選抜ラベルを回帰の応答に使わず、独立実験の植物バイオマスを使用する。

## 正真正銘の実験と外部形質

Anic et al. (2022), [Ecology and Evolution DOI 10.1002/ece3.9157](https://doi.org/10.1002/ece3.9157)、元 [Figshare 24p数値データ DOI 10.6084/m9.figshare.20304039](https://doi.org/10.6084/m9.figshare.20304039)、固定PDF SHA256 \`baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0\`。元真正PDFの前半14ページから**665植物種×原トラフ行、13生存植物種、54真正の物理メソコスム**。3つの光コード（L1/L2/L3）×3つの栄養段階（T1/T2/T3）、9セル各6トラフの原実験を真正検証。13種や665行は独立照明反復ではない。光コードは今回の原データだけでは未照合のため、暗・琥珀色・白LEDと**推測して置換しない**。

植物の湿潤性は独立 [Tichý欧州EIV原資料](https://doi.org/10.1111/jvs.13168) のMOISTを二名法完全一致でのみ結合。真のCSV SHA256 \`ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba\`。別名、集約種、欠測を低湿潤（乾性）に自動割当しない。

元の全54トラフ実験の光コードごとの元平均・SD・SEを利用し、栄養3条件に等重みを付与した各種別の総乾燥バイオマス比（L2/L1、L3/L1、L3/L2）を算出。主仮説の湿潤性との関係は種別Spearman、3コード比較の多重検定は最大絶対ρのMOISTラベル置換（4,999回）を準備。ただし、推定の前に以下の最低共通支持条件を固定：真正MOISTを持つ植物≥10種、そのうち湿潤性>7が≥4種、乾性≤7が≥4種、湿潤性指標レンジ≥3点。共通支持が足りなければ**係数推定の結果には昇格しない**。

## 実際に行った原著PDF＋独立EIV結合結果

[原資料・解析・テスト成功 Actions #44662225191](https://github.com/zuizui0223/alan/actions/runs/44662225191) と [真正数値・種別ソースゲート](../results/all_species_forest/biotic_interactions/anic2022_original_13species_eiv_moisture_lightcode_20261120.json)。

| 実測の共通支持 | 数 |
|---|---:|
| 元PDFに物理トラフ応答がある植物 | **13種** |
| 同じ二名法でMOISTの正数値が得られる植物 | **5種** |
| MOIST>7の真正湿潤種 | **1種** (*Lotus pedunculatus*, MOIST=7.9) |
| MOIST≤7の比較植物 | **4種** |
| 3つの光コード対比をすべて計算できるMOIST測定種 | **4種** |
| 多種MOIST×照明モデレータの採否 | **SOURCE_MOISTURE_MODERATOR_COMMON_SUPPORT_STOP** |

同種内の照明コード別光応答は真の元PDF数値から算出できるが、湿潤植物は実質一種類しか含まない。L3/L1で5種のMOISTとバイオマス比の記述的Spearman ρ=−0.8を得たが、**湿潤種1種類による狭い共通支持・処理色未対応・54トラフ間の種共有共分散未モデル**のため、これを有意な水分嗜好×ALAN交互作用や機構的な結果とは判断しない。光対比方向の代表値を本文で選び直すこともしない。

### なぜEIVの完全一致は5種だけか

[追加原本の分類ラベル監査成功 #44662579384](https://github.com/zuizui0223/alan/actions/runs/44662579384)、[元EIVラベル・欠測理由](../results/all_species_forest/biotic_interactions/anic2022_13plant_unmatched_eiv_taxon_source_labels_20261120.json)。原論文13種のうち8種では真正MOISTの完全一致がない。その内訳に、学名対応があるがMOIST欄自体が欠測の例 (*Dactylis glomerata*, *Holcus lanatus*, *Hypochaeris radicata*, *Trifolium pratense*) と、原EIVが集約名を使う例 (*Achillea millefolium aggr.*、*Lotus corniculatus aggr.*、*Anthoxanthum odoratum aggr.*)、古い同物異名の可能性がある例 (*Agrostis tenuis*) が混在する。例として *Achillea millefolium aggr.* のMOIST 4.1、*Lotus corniculatus aggr.* のMOIST 4.0 は集約群の値であり、元の種の個別実測形質と確認されるまでは**代入しない**。*Anthoxanthum odoratum aggr.* にはMOIST値がない。

これらの分類候補を機械的に採用しても湿潤種を多数増やせる保証はなく、結果が生態学的に肯定になるから採用することは禁止する。最も重要な方法論上の発見は、一般的な欧州草地13種でさえ、**異なる元資料の種×MOIST×照明による本格的な比較の共通支持が不足する**という点である。

## 本文に使える結果と使えない結果

- 元のバイオマス種別光対比（54トラフ、種別平均・ばらつき）自体は実測原資料である。ただし本研究の新しいMOISTとの交互作用係数を保証しない。
- この実験は直接植物体の応答・競争を見ており、*Lu 2026*の5明暗ペア送粉者、*Page 2021*の33研究花粉機能とは別の作用段階である。複数原著が同じ光設備を共有するAnic/Bennie2018/2015は**一つの54トラフ実験プログラム**であり、独立実験3件としてメタ解析に重複計上しない。
- 元のMOIST偏りの主結果は、既存結果として保全できる。今回の機構交互作用を支持できなかったことは、元の偏りが誤りだという検定結果ではない。
- 機構の独立した新しいResultsに昇格するには、水分嗜好が広く異なる植物種を同一照明実験に多数持つか、各研究内で湿潤／乾性種の照明応答を比較できる独立野外実験を5件以上集める必要がある。その際、実測植食・送粉・花量・種子の応答を混ぜずに段階別に解析する。

## 再現

- [元PDF原データの真正復元](../scripts/analysis/alan_anic2022_original_pdf_biomass_54trough_parse_20261119.py)
- [13種×3光コード×元MOISTの原実験検証](../scripts/analysis/alan_anic2022_13species_eiv_moisture_lightcode_source_gate_20261120.py)
- [原本EIVの不一致・集約種と欠測の判定](../scripts/analysis/alan_anic2022_13plant_eiv_name_support_audit_20261120.py)
- [3つのソースゲート合成テスト](../tests/test_alan_anic2022_13species_eiv_moisture_lightcode_20261120.py)
- [GitHub Actions workflow](../.github/workflows/alan-anic13-eiv-moisture-20261120.yml)

**結論：主結果を維持し、確証的な湿潤性×ALANの直接植物機構は source support STOP。現在の資料を超えて肯定的メタ係数を創作しない。**
