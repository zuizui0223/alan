# 夜間送粉と自動自殖を同じ種で照合する

[PR 84](https://github.com/zuizui0223/alan/pull/84)の自動自殖と[PR 85](https://github.com/zuizui0223/alan/pull/85)の昼夜排除実験を種名で結んだ。PR 85 の先端 8c0c610ddd3e00879fa20b2220425184d4c0da45 の元モデルと公開実行成果物を使い、52種中39種の推定と平均群間差 −1.646 pp を再現した。10島以上の22種のうち、夜優位かつ自動自殖なしは1種、自動自殖ありかつ夜優位は0種。したがってこの支持範囲で両形質の交互作用は識別できない。

- [全実験134種](species.csv)。夜優位とは文献の夜のみ／（昼のみ＋夜のみ）繁殖成績の種別集約値が0.5超という意味。実測の自然下の夜間寄与率、効果の有意性、蛾の独占を表さない。欠測は不在を表さない。自動自殖は送粉者なしでの繁殖能力を表し、自家和合性や遺伝的自殖率とは異なる。
- [9候補と元の78地域](candidates.csv)。同じ5年・10,698地点年に戻ると9種すべてで島の記録0。島・本土共通544種のうち実験があるのは13種。夜優位の3種では S. latifolia と S. nutans は両側負、T. pratense は両側正。3種とも島−本土差は負。少数地点の写真構成比傾きである。
- [PR 85 の島52種と推定値](pr85_island_species.csv)。島年内の先行光勾配に対する写真記録確率の変化（pp）。自動自殖なし・夜優位は5種、このうち4種で推定できた。G. conopsea −1.351（9島）、S. stellata −0.941（2島）、P. persica −1.870（28島）、I. vera +4.098（4島）。I. ingoides は支持不足。正負とも残す。種ごとの区間や新たな検定は行っていない。島数は種ごとに情報を与えた島数。元の島−本土モデルの現年VIIRS・群集構成比とは尺度も地域も異なる。

G. conopsea は自家和合でも媒介者が必要で、スウェーデン4集団とノルウェー2集団では昼と夜の重要度が逆転する。追うべき尺度は種内の地域集団である。[地域差と花香](https://doi.org/10.1093/aob/mcx203)。
S. latifolia は雌雄異株だが直接形質台帳では欠測のため9種の外に保持した。51島の写真傾き −0.376。独立した野外照明試験では同じ蛾が送粉し幼虫が種子を食べ、照明で被害果実の割合が変わっても推定種子数の差は不確かだった。[Boom et al.](https://doi.org/10.1038/s41598-020-75471-1)。
固定78地域の既存監査では野生成虫の攻撃・捕獲を努力量つきで島と本土の明暗で比較できる地域は0。9候補はいずれも既存の両側開花11種に含まれない。地点年での機序は未検証である。

## 再現

PR 85 の上記版の scripts を展開した作業場所を PR85 とする。元ワークフローの署名済み入力を使い、PR85 の下に pcnl、平置きした groups、annual/yearYYYY、climate、ghsl/ghsl2015_2020 を置く。以下は順に134種、固定地点照合、PR 85 の推定を再作成する。生の座標は出力しない。

    python scripts/analysis/alan_joint_night_autonomy.py PR85 AUTHOR_DATA DIRECT_TRAIT_LEDGER_GZ
    python scripts/analysis/alan_joint_candidate_sites.py FROZEN_V011_PICKLE_DIR STRICT_SITE_LIGHT_DIR
    python scripts/analysis/alan_joint_pr85_islands.py PR85 SIGNED_GLOBI_GZ

[2023年の全体平均の比較対象](https://doi.org/10.1002/ece3.9990)。主要結果と64種の副次結果は既存の定義を維持した。
