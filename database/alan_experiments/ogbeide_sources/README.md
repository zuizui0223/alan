# Ogbeide2026 aquatic–terrestrial food-web experiment

Primary article: https://doi.org/10.1111/1365-2435.70335 . Original data: https://doi.org/10.5281/zenodo.14516568 .

## Design appraisal

The six-week Landau experiment crossed light and crayfish in 16 stream–bank mesocosms, four per cell. Emergence samples are nested in flumes. Consumer isotope samples were taken later than prey samples. Companion Manfrin2025 may share the programme. Isotope mixing proportions are model estimates, not attack rates. Published uncertainty for ALAN-related changes in spider aquatic dietary contributions spans zero. A change in isotopic niche width is not evidence of increased total consumption, moth predation or reduced pollination. Retain dietary composition, emergence and predator abundance as separate outcomes. Source: primary Methods2.1–2.8, Results3.1; article model not re-fitted here.

## Original-data acquisition and audit

Six individually deposited CSV files and the 322,650,517-byte code/data archive were retrieved. Zenodo MD5 values were verified. Relevant archive members were checked through ZIP CRC and compared byte-for-byte with the initially read entries. Files remain in the external cache. `download_manifest.csv` gives hashes and URLs. Individual and archived isotope CSVs have different byte sizes; archived files contain greater numerical precision. Do not interchange them silently.

Our audit of the standalone consumer tables found 48 spider measurements distributed across 16 flumes (three per flume) and 48 crayfish measurements across eight flumes (six per flume). These are not 48 lighting replicates. No isotope mixing-model replay or diet effect-size extraction has yet been performed.

## Recalculated emergence snapshot

`emergence.rds` contains 192 records: two sampling weeks × 16 flumes × three traps × two collection labels. Both labels carry four trap-days and unit area in the supplied table. This recorded exposure is retained; calendar coverage and any overlap need clarification. The stored rate equals count / trap-days / area in every row. Each flume has six records per wave with equal effort, so its mean rate also equals total count / total recorded area-time.

The latest available wave is Week4, while the experiment continued to Week6. `extract.py` exports flume rates and two Week4 ALAN/dark comparisons, keeping crayfish status fixed:

| Crayfish | Lit mean | Dark mean | Ratio | Approximate 95% interval |
|---|---:|---:|---:|---:|
| Absent | 22.7813 | 33.1979 | 0.686 | 0.261–1.803 |
| Present | 18.2708 | 12.8854 | 1.418 | 0.424–4.740 |

Units are individuals per square metre per recorded trap-day. Each arm has four flumes. Intervals use independent-flume delta variance and a Welch-type t approximation, not the article's temporal model. The difference between light lnRR with and without crayfish is +0.72576 (sampling variance0.39832); these point estimates do not establish an interaction. Both light-effect intervals include one.

The combined ALAN+crayfish arm versus untreated control is a different estimand from ALAN+crayfish versus crayfish alone. Do not attribute the whole combined-stressor contrast to light. The two extracted comparisons stay outside the main final-time catalogue because Week6 data are absent here. They do not alter the wet-affinity macro result.

## Replay

From repository root, export the original object using R `readRDS` to `selected/emergence_export.csv`, then run `python database/alan_experiments/ogbeide_sources/extract.py`. Input SHA256 values are in `source_manifest.csv`. The script checks exposure calculations, unique observation keys, treatment identity, replicate counts and source endpoint. No zero correction or model-based standard errors were substituted.
