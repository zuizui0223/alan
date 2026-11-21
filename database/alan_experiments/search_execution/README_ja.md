# 契約検索の実行・受け渡し

準備日：2026-11-21。利用者は両データベースを利用可能と回答済み。
こちらのブラウザ操作機能が初期化に失敗したため、契約セッションの確認・検索・エクスポートは未実施。
検索式を作成した日を検索実施日として記録しない。既存のPRISMAログの未実施状態を維持する。

## 実行

1. Web of Scienceで **Core Collection** を選び、Advanced Searchへ
   [WoS検索式](wos_advanced_query.txt)をそのまま貼り付ける。
2. ScopusのAdvanced Searchへ[Scopus検索式](scopus_advanced_query.txt)を貼り付ける。
3. 期間は利用可能な全期間。言語・分野・文献種別・Open Access・湿地などで追加絞り込みをしない。
   実際の検索日、検索結果総数、WoSで利用できた索引／年代範囲を記録する。
4. 全結果を出力する。一度に全件を出せない場合は、画面の出力上限に従って連続した範囲に分け、各範囲を記録する。

エラーになった検索式は黙って短縮しない。エラーと修正版を記録し、元の概念が保たれるか確認する。
実行確認済みの式ではなく、公式のフィールド指定に合わせた実行候補である。

## 受け取るファイル

- **WoS：Plain Text（テキスト）**。レコード内容はFull Record、可能ならFull Record and Cited References。
  題名、著者、抄録、年、誌名、DOI、WoS識別子を含める。
- **Scopus：CSV**。Citation information、Bibliographical information、Abstract and Keywordsを選ぶ。
  DOI・EID・題名・著者・年・誌名・抄録・キーワードを含める。利用可能なら引用文献も含める。
- 各検索の実際の件数・日付・設定を[記録表](search_execution_template.csv)に残す。
  この表の空欄は未実施・未取得であり、ゼロ件ではない。

元ファイルは編集せず渡す。受領後にファイルのハッシュ、総行数、識別子重複、抄録の欠落、
全件出力との一致を確認する。元エクスポートはローカル資料として保管し、公開リポには自動追加しない。
現在の553報告台帳へDOI・索引IDで照合し、曖昧な題名一致だけで統合しない。
新規検索の件数が確定してからPRISMAログと選別対象を更新する。

## 仕様の根拠と検索対象の違い

WoSのTSはTopic、ScopusのTITLE-ABS-KEYは題名・抄録・キーワード検索。
両者の索引・フィールド範囲は完全には同じでなく、件数を一致させるために式を調整しない。
既存プロトコルの3概念群と語は維持した。

- [WoS公式：Advanced Search field tags](https://webofscience.zendesk.com/hc/en-us/articles/26916347018257-Web-of-Science-Core-Collection-Advanced-Search-Field-Tags)
- [WoS公式：検索結果の出力](https://webofscience.zendesk.com/hc/en-us/articles/20135124532113-Search-Results)
- [Scopus公式：Advanced Search](https://www.elsevier.support/scopus/answer/how-can-i-best-use-the-advanced-search)
- [Scopus公式：エクスポート項目](https://www.elsevier.support/scopus/answer/how-do-i-export-documents-from-scopus)

これらは操作・形式の参照先であり、実際の検索を行った証拠ではない。
