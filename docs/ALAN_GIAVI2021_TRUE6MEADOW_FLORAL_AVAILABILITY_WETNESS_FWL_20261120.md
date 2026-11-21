# ALAN湿潤植物の機構Results：真正Giavi2021 6草地の花量・昼間訪花・湿潤性の関係

**重要な後続監査（2026-11-20）：この文書のPLANTATT花被度β_F=−0.769・名目符号反転p=.03125は、湿潤種固有のALAN効果の証拠には使えません。** [原本の同一草地×同一調査回に絞った反証](ALAN_GIAVI2021_SAME_MEADOW_WETDRY_SELECTION_FALSIFICATION_20261120.md)では、真の湿潤種と比較種の同時比較が **PLANTATT 3調査回／2草地、Döring 7調査回／2草地** に限られること、湿潤種の同一回明暗両側への記録割合が **8/27=29.6%** と比較群の **69/130=53.1%** より低いこと、湿潤種2種を完全除外してもF≤7内で負の花被度差勾配が残ることを実測確認しました。**この文書の回帰係数は探索的な異種間・草地間関連であり、光による湿潤植物の開花減少機構は識別できていません。** 測定対象は花数でなく調査線に沿った花の被覆面積プロキシです。

2026-11-20 | Source-verified independent raw-data reanalysis | Draft only | No wet-plant causal mediation claim

## 研究の接続

Mainの最初の発見は、長期の増光地点で記録率低下候補に選ばれた植物群の湿潤性が高いこと（環境調整後 +0.534）。その原因が「昆虫の訪花選好の変化」なのか「花の利用可能な量の変化」なのかは未確定。これは主結果の種群構成と別の独立した同時期の野外照明実験で検証する機構的候補の一つであり、主結果の因果媒介効果を推定するものではない。

[Giavi et al. 2021 Nat Commun 12:1690](https://doi.org/10.1038/s41467-021-22011-8) original published [Zenodo RData](https://doi.org/10.5281/zenodo.4540407) MD5 f60b21d82c163cbb53520fda0bc61985. Exact original RData 1,884 observed treatment×site×sampling×plant×insect-order rows, 6 physically paired lit/dark meadows, 37 paired survey events, 21 plant label, three source guilds (Diptera, Hymenoptera, Coleoptera). 1,441 original observed zero guild visits, 0 missing source visits. Same plant in same true meadow pair and survey with both C/T original observations: 193 plant×site×sampling pairs. We do NOT assign 0 to unreported plants.

External plant moisture affinity: RBG Kew BIFloraExplorer original data-raw/BI_MAIN.csv, source Git SHA1 e6a1da929af81f370d2a1b79177b489e3e573a06, separate British PLANTATT F and Döring F scores. Independent original taxon names recognized exactly, no automatic assignment of genus-only/aggregate/NA.

The original authors had already controlled flower abundance in their insect-group mixed-effects model. Our *new* test asks whether the **MOIST indicator of a species** changes its daytime visitor response and how much of the association is statistically aligned with *concurrent flower-cover differences*. This is an exploratory reanalysis, not a re-publication of the original mean treatment effect.

## Original-data new methods

- For every original site×sampling event×plant×arm, the three insect guild-specific source rows are required to agree exactly on measured original flower proxy plab; visits across guilds are summed but flower cover counted **once**.
- Source Δ log1p(day visits) = log(1+observed lit visits) − log(1+observed dark visits); source Δ log1p(flower proxy) similarly. It is not a fruit/seed or true stigma-compatible pollen delivery measure.
- Fit ordinary descriptive sampling-event fixed-effect coefficients on source-matched plant observations (6 site pairs, 37 sampling events); compare F slope before vs after **concurrent** flower-cover difference as covariate. Separately fit flower-cover difference as an outcome. Three original insect guilds are supplementary different endpoints.
- The algebraic Frisch-Waugh-Lovell identity from identical source rows is verified: initial β_F = β_F adjusted for flower + γ_flower→visitation × β_F→flower. This is an exact **association partition**, NOT a path model or an estimated mediator attributable to ALAN.
- Descriptive exact signflip reassigns 6 original meadow C/T labels in whole-site units (64 possibilities), not independently for the 193 plant observations or 37 visits. This is **not an identified randomized-treatment test**; the source paper does not provide a verified randomized site-assignment ledger.
- Source **common-support stop**: we require ≥4 distinct wet plants (F>7), ≥4 drier plants, and ≥3 true meadow pairs with ≥2 wet and ≥2 drier scored plants within the same original meadow. The actual source fails this for both F indicators.

### Source-verified main descriptive results

[CI #44669315284 SUCCESS (final source and whole-plant leave-one-out)](https://github.com/zuizui0223/alan/actions/runs/44669315284). [Verbatim full aggregate source output](../results/all_species_forest/biotic_interactions/giavi2021_original6meadow_flower_FWL_wetness_20261120.json).

| Original separate trait | PLANTATT F | Döring F |
|---|---:|---:|
| Source scored source species | 7 | 8 |
| True matched C/T plant×sampling pairs | 53 | 77 |
| Source wet (F>7) species | 2 | 2 |
| Distinct meadow sites with ANY wet+dry | 2 | 2 |
| Distinct meadow sites with TWO wet+TWO drier | 1 | 2 |
| β_F for lit-minus-dark log visits, without flowers | **−0.302** | **−0.412** |
| β_F for lit-minus-dark log visits, with flower difference | **+0.098** | **−0.016** |
| β_F for lit-minus-dark log flower cover itself | **−0.769** | **−0.646** |
| Algebraic component associated with concurrent flower differences | −0.400 | −0.396 |
| Whole-meadow label signflip nominal absolute tail for adjusted visit slope | 0.344 | 0.969 |
| Whole-meadow label signflip nominal absolute tail for flower-F slope | **0.03125** | 0.250 |

The apparent negative wetness–visitation slope becomes ~0 for the all-source-scored cohort with a contemporaneous flower-cover covariate. The flower-abundance difference itself shows a negative moisture association, with **indicator sensitivity**: PLANTATT's conditional meadow signflip p=0.03125 versus independent Döring F p=0.250. No multiplicity or source group selection correction is supplied for these exploratory contrasts; the threshold alone is NOT independent confirmation of a wetness-specific mechanism.

### Repeated-meadow robustness and within-meadow support

Restrict to plant species documented with both C/T source arms in two or more *different* meadow pairs:
- PLANTATT F: **4 species / 42 plant×event pairs / six source sites**; visit β_F **−0.409 before and −0.401 after flower adjustment**, whole-site signflip for adjusted visit p=0.1563; flower β_F **−0.817**, whole-site signflip p=0.03125.
- Döring F: **5 species / 66 pairs / six sites**; visit β_F **−0.385 before and −0.140 after adjustment**, whole-site signflip for adjusted visit p=0.4688; flower β_F **−0.573**, whole-site signflip p=0.25.

This restriction is a critical counterexample to claiming flower resources explain *all* of the originally observed wetness-visitation association: on the tiny repeated-plant subset, PLANTATT visits remain negative after floral adjustment. The original source still has only 2 true meadow sites co-observing wet and drier plants and cannot identify a general wetness moderator.

With only same-meadow source sites that have both wet and drier plants (two physical sites total):
- PLANTATT: 5 plants / 14 plant×sampling pair rows, visit β_F −0.211 before, +0.009 after floral adjustment; corresponding flower β_F −0.644. Two-site signflip can only be coarse (4 original assignments), p=1 for adjusted visits.
- Döring: 6 plants / 22 rows, visit β_F −0.098 before, +0.169 after floral adjustment; flower β_F −0.479; two-site flip adjusted visits p=1.

Thus the descriptive association may be driven by **flower resource availability** or correlated site/plant composition; the design cannot distinguish those paths from a treatment-induced change in flower production/phenology. In the true matched source, the main flower-abundance correlation is not immune to cross-site confounding or different phenological survey coverage.

### 最終感度分析：元の植物1種を丸ごと除く

[真正原本からの最終source CI #44669315284（成功）](https://github.com/zuizui0223/alan/actions/runs/44669315284) では、植物種ごとにすべての調査日・草地記録を除去してから係数を再推定した。花量差への湿潤性の負方向はどの植物を除いても維持される。

| 元形質 | 全植物の花量差 β_F | 1植物種除外後の β_F 範囲 | 除外した種の数 |
|---|---:|---:|---:|
| PLANTATT F | −0.769 | **−0.933〜−0.641** | 7 |
| Döring F | −0.646 | **−0.847〜−0.506** | 8 |

したがって花量との対応が**単独の1植物種のみ**で生じたとは考えにくい。ただし共通支持のある湿潤種は元々2種・真の湿乾共観測は2草地なので、「全湿潤植物に普遍的な効果」とは認めない。一方、**花量を調整した後の訪花数 β_F は、1種を除くと正負が入れ替わる**（PLANTATT −0.165〜+0.612、Döring −0.097〜+0.139）。この不安定性こそ、訪花者固有の湿潤性感受性をResultsで主張しない理由である。

## How to place it in manuscript

Keep Results 1 = robust long-term observation of wet-plant bias. A possible independent **exploratory Results** panel may state:
"In a separately manipulated six-paired-meadow system, apparent negative relationships between species moisture affinity and the illuminated-minus-dark diurnal visitation response were substantially attenuated when measured floral availability was included. In parallel, the moisture gradient was negatively associated with the relative floral-cover difference, suggesting that plant resource availability and insect visitor responses should not be conflated. However, only two meadow pairs provided simultaneous wet and drier plant comparisons, and the wetness moderator was not identifiable as a generalized light-specific causal pathway."

Do NOT call this a confirmed ALAN effect on plant flowering or a measured mediation percentage. In Discussion, retain pollen-transport efficiency, nocturnal visitor replacement, grazer pressure and wetland drainage as unresolved alternatives.

## Reproduction

- [Source-authentic visit/flower extraction and meadow-preserving FWL](../scripts/analysis/alan_giavi2021_true6meadow_floweradjusted_moisture_20261120.py)
- [Synthetic source cluster/fixed-effect and decomposition tests](../tests/test_alan_giavi2021_true6meadow_floweradjusted_20261120.py)
- [GitHub Action source-replay workflow](../.github/workflows/alan-giavi-true6meadow-floweradjusted-20261120.yml)
- [Complete source aggregate report](../results/all_species_forest/biotic_interactions/giavi2021_original6meadow_flower_FWL_wetness_20261120.json)

Scientific next data need: independent multispecies lit/dark experiments providing within-meadow wet/drier plant contrasts across many independent meadow pairs, **flower cover and compatible stigma pollen/seed endpoints**, and true pretreatment floral abundance. Such original data can separate loss of flowers from visitor redistribution; external broad wetness macro link still needs independent validation.
