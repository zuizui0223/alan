Current update: 146 retained reports reviewed,42 pending; nine outside-search reports. The final ALANizer biological paper's indexed Methods/Results are now read. Its216samples do not match268near/far values in the public table; exact inclusion mapping remains unresolved. Final methods-paper version still unread. No effects added.

Previous update: 146 retained reports reviewed,42 pending. ALANizer v2methods/validation and linked134-row biomass source audited;8records differ from planned allocation but all match corrected metadata. Final methods and biological articles identified as two outside-search leads (nine total); neither final Methods/Results read. No effects added. Shredder2021/2022 remain access holds.

Previous update: 145 retained reports reviewed,43 pending. Odonate polarization preprint and final Methods/Results reviewed; both excluded from ALAN manipulation synthesis because observations used natural sunlight at11-14h. A third within-search version pair yields514 publication families; retained families186 and reports188 unchanged. Final species models weaken several preprint preference claims. No pooled effects added; final eligibility remains unfinished.

Previous update: 144 retained reports reviewed,44 pending. Daufel original workbook recovered and independently checked; two conditional oviposition estimates remain outside the catalogue because daily completeness and single-egg total differ from the paper. Existing pooled results unchanged.

Current update: 143 retained reports reviewed, 45 pending. Castrop plant preprint reviewed; ten author-reported percentage contrasts lack extracted variance and remain outside the finite-effect catalogue. Raw access and shared lighting-row support unresolved.

Current update: 142 retained reports reviewed, 46 pending. Two verified within-search version pairs yield515 publication families (186 retained); reports remain517 (188 retained). Moth radar versions are one experimental programme, with different endpoint denominators. No pooled effect added.

Current update: 141 retained reports have primary Methods/Results reviewed; 47 remain pending. Biofilm v2 is reviewed but subsample variances are not eligible for pooling. Final biofilm article remains unread; outside-search report count is seven. No new finite effects.

Current update: 140 retained reports have primary Methods/Results reviewed; 48 remain pending. Wang grassland preprint and published version reviewed separately; final DOI is an additional outside-search report (six total), not an additional experiment. Raw numerical extraction remains pending.

Current update: 139 retained reports have primary Methods/Results reviewed; 49 remain pending. Amphibian2018 final-paper review completed; two freshwater predator preprints remain access holds with probable overlap, not merged. Final eligibility and PRISMA remain unfinished.

# 検索記録と本文確認の照合（2026-11-20）

**最新更新（Abonyo/Oduor2024確認後）**：対象188報告＝本文確認138＋未確認50。AGR:IND608156393を巻号・ページ・著者で同定し、方法・結果を確認した。研究数517・出版系列516・検索レコード521は不変。下記の137／51は直前の履歴。

**最新再集計**：521検索レコード、517個別報告、確認済み版対応後516出版系列。全文確認対象188報告＝本文確認根拠あり137＋未確認51。最終採択・独立実験数ではない。AGR:IND608754560と既存MED:39692005を一次本文・AGR著者順・巻号・年で同定したため、未確認数が1減った。新たに1報読んだという加算ではない。Weschke2024では公開原ファイル23件を取得し、最終20地点の食性群別記録を再集計した。検索外の学位論文は別の`grey_literature_screening.csv`に保持し、既存DOI照合表の検索外5報へ自動加算していない。以下は初回と途中の照合履歴であり、現状数には最新JSONを用いる。

補助検索の521レコードは、現時点で識別子により区別した519論文グループに対応する。タイトルだけの重複候補やプレプリントと最終論文の関係は未解決なので、519を最終的な独立論文数とはしない。実験系列数とも異なる。

一次選別で全文確認へ残した192レコードは、確認済み重複をまとめると190論文グループになる。そのうち98件には、本文レビューまたは原データ抽出を伴う本文確認の根拠がある。残り92件は本文確認の根拠が不足する。アクセス保留、部分的な方法確認、書誌のみの状態を区別したまま`pending_primary_review.csv`へ残した。98件は採択数ではなく、確定メタ解析標本数でもない。

従来のバッチ台帳94件（本文確認90・アクセス保留4）に、別途記録した追加監査や抽出を照合した。ミミズ実験DOI10.1186/s12862-024-02200-xはバッチ7と以前の土壌監査の両方で確認されていた。したがってバッチ7の4件を、新たに読んだ独立4論文として加算しない。履歴は消さず、同じ論文として扱う。

AGR:IND609539107はDOIもPMCIDも欠けていたが、EuropePMCの原書誌と一次論文のタイトル、著者順、雑誌、巻号、年を照合し、DOI10.1111/gcb.71029・PMC13421007と同定した。既存のMEDレコードと同じ論文で、既にバッチ6で本文確認済みだった。`identity_adjudications.csv`に判断と一次出典を記録し、原AGRレスポンスを保存した。自動タイトル一致だけで統合していない。

検索外の追加出典1論文も別表に保持した。一次検索と引用追跡等から得た論文を混ぜて検索ヒット数を増やさない。Parkinson2023とDeitsch2023を研究登録表へ追記し、抽出済み台帳の全9論文学名ラベルが研究登録表に存在することを検証した。

検証：521レコードの保存、190=98+92の分割、同定したAGRとMEDの同一論文への結合、当該論文の本文確認への接続を確認した。これは記録照合の検証であり、独立した第二査読者による採択判断ではない。WoS/Scopus検索、残る本文確認、適格性確定、効果量の第二確認、段階別メタ解析は未完了。

再現：`python database/alan_experiments/reconcile_search_identity.py`、続いて`python database/alan_experiments/reconcile_review_coverage.py`。出力は`report_review_coverage.csv`、`pending_primary_review.csv`、`supplementary_report_sources.csv`と各集計JSON。

追加更新：Ogbeide2026の一次本文と公開原データを確認し、AGR:IND609471227を著者順・巻号・ページを含む元書誌と照合した。未確認93件。バッチ台帳94件とは集計範囲が異なる。

追加更新：Grognuz2026の本文と公開第1版を確認。第2版との対応と実験系列の重なりは保留。現在の本文未確認は92件（採択保留件数とは異なる）。
