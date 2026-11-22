# Sullivan et al. (2019) 2013湿地実験：真正原データの新規再解析（Draft 2026-11-22）

**独立した新しい実験プログラム1つの候補。統合メタへの採択・湿潤性メタ回帰は未実施。**

**原著**: Sullivan, Hossler & Meyer (2019), *Ecological Applications* 29:e01821, DOI [10.1002/eap.1821](https://doi.org/10.1002/eap.1821)。**原本**: Dryad DOI [10.5061/dryad.gb9534p](https://doi.org/10.5061/dryad.gb9534p)。公式Dryad REST v2のversion 8787の公開ZIPから真正原ファイルを抽出し、**観察研究2010–2011年の別ファイルを除外**した。2013実験の元TSV SHA256 `385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484`。Dryadがその場で作る外側のZIPのバイトハッシュはダウンロードごとに変わり得るため、**検証対象はZIPそのものではなく内部原TSVのSHA256**とし、外側ZIPのハッシュは実行別の受領記録に残す。元TSVの原行・希少地点情報は本PRに再配布せず、取得スクリプトとハッシュ・派生統計を保存。

## 真の単位と抽出

- 2013年の**1,161原分類群・採集記録／13列**。単位を厳格に復元すると**15の物理ステーション**＝3 habitat/location `Group` ×5照明処理（各ブロックに暗区DRK・LEDI・LEDH・HPSI・HPSHが1ステーションずつ）。各処理の物理割付反復は**3ステーション**で、分類群行や採集日が独立した照明処理ではない。
- 各処理と暗区が同一`Group`×同一`Date`でそろう日時は、水生羽化昆虫トラップ（`AquaticEmergent`）、陸生節足動物（`TerrestrialArthropod`）、クモの網数（`WebCount`）で各ブロック**4採集日**。原表の観測行をイベント単位で合算し、原採集日欠損を「捕獲0」に置換しない。葉リター（`LeafFall`）はCount非該当の別応答で、昆虫捕獲へ混ぜない。
- 原表の個体捕獲数（計数表の記述和、**独立標本数ではない**）は、水生羽化昆虫 DRK=407、LEDI=1212、LEDH=1362、HPSI=1565、HPSH=826、陸生節足動物 DRK=161、LEDI=205、LEDH=179、HPSI=125、HPSH=169、網数 DRK=226、LEDI=66、LEDH=132、HPSI=58、HPSH=183。比の主要定義はこの集計和の単純比ではなく、**3独立Groupごとに4日分を集約した照明／暗区比を対数化し、その3つを平均したもの**。

## 探索的な明暗応答

| 原応答 | LEDI RR | LEDH RR | HPSI RR | HPSH RR |
|---|---:|---:|---:|---:|
| 水生羽化昆虫捕獲 | 2.49 | 3.04 | 3.23 | 1.83 |
| 陸生節足動物捕獲 | 1.25 | 1.11 | 0.79 | 1.08 |
| クモの網数 | 0.21 | 0.57 | 0.23 | 0.57 |

各対比のSE、**原Group間の4×4共分散**、自由度2のt係数に基づく**記述的**95%区間は [experimental_endpoint_lnRR_candidates.csv](experimental_endpoint_lnRR_candidates.csv) と [真正再実行JSON](dryad_provenance_and_frozen_block_results.json) を参照。例えば、水生羽化昆虫の**LEDH/暗区**はRR **3.04**（記述的95%区間 **1.05–8.82**）、網数の**HPSI/暗区**はRR **0.23**（**0.073–0.728**）。少数3ブロック・多重探索であり、事前確認的有意差の主張に使わない。

**同一実験内の相対反応**も、同じGroupで得た差（`lnRR_aquatic − lnRR_terrestrial`など）から計算した。水生捕獲の対陸生捕獲の相対比はLEDI **1.99**、LEDH **2.73**、HPSI **4.08**、HPSH **1.69**。LEDHの記述的区間は **1.17–6.33**。水生捕獲の対網数はLEDI **11.93**、LEDH **5.36**、HPSI **14.05**、HPSH **3.20**。全区間とブロック別原値は [within_experiment_cross_endpoint_ratios.csv](within_experiment_cross_endpoint_ratios.csv) へ。これらは交差段階の*観測数反応差*であり、**捕食率や捕食による種子生産の媒介推定ではない**。

## 採択・因果の境界

1. **実験プログラム単位の新規候補は1つ**。15地点・3生態段階・4光色を複数の独立研究としてメタ統合しない。原著の照明割当て無作為化・光漏れ・空間交絡・標準化採集努力の完全な人手監査が未完なので、共通効果量カタログ（108行）やmeta kに**未追加**。
2. 同一湿地内で3 habitat/location groupを比較したもので、**湿潤地と乾燥地の対照**ではない。植物種EIV、水面距離×照明、植物繁殖はこの原表にない。ALANによる湿潤植物の記録率低下の原因は同定できない。
3. 網数の低下はクモ個体数や実捕食率の低下と同義ではない。水生羽化昆虫の捕獲増加は実出生・個体群増加の証明ではない。光色・強度・複数反応を事後的に探索したので、多重比較補正後の確認的効果とは呼ばない。
4. 同一研究の原標本で3応答が同時に測定されても、機構鎖の因果媒介係数や植物への効果は推定できない。他のParkinson 2020/Manfrin 2025と相反する捕食者応答は、**文脈依存の検討対象**であって、どちらかの研究を排除する理由ではない。

GitHub Actionsによる真正ソース再取得と10 unittestは [ワークフロー](https://github.com/zuizui0223/alan/blob/analysis/alan-arthropod-search-expansion-20261122/.github/workflows/alan-sullivan2019-dryad-source-gate-20261122.yml) による。コード `scripts/acquisition/alan_sullivan2019_dryad_experiment_source_gate.py` は再実行で原ハッシュ、群・日・照明処理の比、共分散を復元する。派生CSVは原表の個別観測を再配布しない。
