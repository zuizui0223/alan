# ALAN — global island alternative-pollinator evidence, and real four-island predator/honeybee floral network

**2026-11-14 original-data report.** This source branch is Draft, stacked on [PR #70](https://github.com/zuizui0223/alan/pull/70), with parent source findings unchanged. The key ecological question is *when* island plants with historically reconfigured nocturnal and/or day visiting channels will be vulnerable to nightlight because functional alternatives, plant phenology, selfing or predator-mediated services do not buffer effective pollen delivery. **This original analysis provides measured island pollinator-channel baseline but no light response.** No raw original sensitive plant×site×visitor observations are published in this repository.

## First: global island/mainland lizard alternative vector evidence tiers

[Justicia Correcher et al. (2023), Ecography DOI 10.1111/ecog.06425](https://doi.org/10.1111/ecog.06425) compiled **452 reported plant–lizard flower-use relationships from 162 primary publications**, involving **173 lizard taxa, 311 plants**, with observations on **91 islands (41 archipelagos)** and in 31 mainland countries. Approximately **79% of all positive relations** were reported on islands. Of 452 published flower-use records, only **90 (20%)** had experimentally confirmed pollination effectiveness, and only **three** plant species had established obligate lizard pollination. These are the original review authors' findings, not independent reanalysis or biological absence estimates. Each reported record is not an independent island, predator encounter, and plant reproductive experiment. The 79% observed island share cannot be converted into a baseline prevalence rate without comparable plant/lizard search effort between mainland and island.

**Raw-source access gate actually executed:** [Actions #43575594251 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43575594251). The authors' 129.5 kB original supplementary interaction workbook, ecog12967-sup-0001-AppendixS1.xls, was located, but all four official Wiley original-source URLs returned **HTTP 403**, so we **did not re-fit** the authors' original 452 interactions or claim new island/mainland efficacy regression. The separate evolution paper [10.1016/j.ympev.2023.107839](https://doi.org/10.1016/j.ympev.2023.107839) has public Figshare [23460260](https://doi.org/10.6084/m9.figshare.23460260), and the original API was successfully retrieved, but that exact archive includes phylogenetic trees/sequence files (3 files) and **not the needed 452 interacting-plant trait rows**. No unrecorded species are coded true biological pollinator absence.

## Second: independent true original Ogasawara four-island floral partner data

Quitián et al. (2026), [Conservation Biology DOI 10.1111/cobi.70304](https://doi.org/10.1111/cobi.70304), original CC BY 4.0 published [Zenodo DOI 10.5281/zenodo.19221853](https://zenodo.org/records/19221853), data file Data_Quitianetal2026_B.xlsx. Original file author-deposited MD5 **68f295b91605ea82bce0c05bffc7f884**, independently verified, with raw original SHA256 **4e0c4a293e2a522e67bc1f87f567249d45b42000b78eaef9d355e8658b1987a6**. Original file exactly 129,502 bytes, one sheet (Data interactions), **2,745** original rows, 11 explicit source fields including island, invasion context, forest condition, anole presence, May/July/September 2021, plant, visitor taxonomic group, visitor identity and count N.Int.

[Original workbook integrity/source schema Actions #43575803117 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43575803117); [full source network and Anijima matched ecological comparison Actions #43576471698 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43576471698), [aggregate-only source artifact #14106935254](https://github.com/zuizui0223/alan/actions/runs/43576471698/artifacts/14106935254).

### A critical original data denominators audit

The author-deposited raw workbook itself includes **5,703** N.Int legitimate/recorded flower-use events, **61** originally labelled plant names and **143** originally labelled visitors/morphotypes, whereas the paper reports **4,174** eligible main-network visits, **55** plants and **115** visitors. Do not silently assert these are exactly equal.

Actual original-family count audit found **1,421 visits from Formicidae ants, 100 from Collembola and 8 from Acari**. Removing these groups gives exactly **4,174** visits (=5,703−1,421−100−8). Dropping a further **25 original rows with N.Int=0** produces **2,063 positive main-network source rows and exactly 55 plant labels**. The authors explicitly state ants were removed from the main networks; the treatment of springtails/mites has not been fully checked from the authors' methods, so the full filter is a **source reconstruction diagnostic**, not an independently pre-registered author preprocessing rule. In the resulting original data there remain **116** visitor name/morphotype labels versus **115** reported in the article, a **one-label unresolved difference**; do not silently delete a pollinator to manufacture 115.

The raw author source encodes the honeybee as **Apis_mellifera** with an underscore, not the Latin space form. The exact label is verified; naïve name-matching would have incorrectly made its proportion zero.

### True source-tested island honeybee floral resource-use gradient

Actual original paper-main-denominator 4,174 daytime flower visitor events (08:30–15:30 only) in 2021, grouped by physical island, original honeybee named positive records, and three seasonal source surveys:

| Physical source island | Source counted flower visits | Apis mellifera visits | Honeybee fraction of source visits | Other measured visitor visits, not necessarily native |
|---|---:|---:|---:|---:|
| Chichijima / 父島 | **1,776** | **1,488** | **83.78%** | 288 |
| Hahajima / 母島 | **631** | **215** | **34.07%** | 416 |
| Anijima / 兄島 | **871** | **95** | **10.91%** | 776 |
| Ototojima / 弟島 | **896** | **4** | **0.45%** | 892 |
| **All 4 physical islands** | **4,174** | **1,802** | **43.17%** | 2,372 |

This huge source-documented **daytime floral visit channel** difference is directly relevant to island pollination-system reorganization. **It is not a confirmed replacement of historic moths**, not directly comparable between islands at fixed taxon composition, and not *effective conspecific pollen deposition*: in this study legitimate visitors touched reproductive organs but single-visit effectiveness and resulting seeds were not measured for every visitor. All non-honeybee positive records are NOT presumed native.

### Anole–pollinator associations are confounded by original island and site

The author original categorical source has **7 invasion contexts** across 4 physical islands (father island disturbed and natural both anole-present; mother island similarly anole-present; Ototojima no anoles; Anijima has one original anole-present and one anole-absent **both natural forest**).

The Anijima source natural forest context-specific result is:

| Anijima source survey context | Actual daytime positive flower visits | Apis source visits | Apis proportion | Mean original plant×season partner diversity, eH |
|---|---:|---:|---:|---:|
| Anole absent, ANI_A | **616** | **45** | **7.31%** | **2.719** |
| Anole present, ANI_P | **255** | **50** | **19.61%** | **2.232** |

The diversity eH is **exp(−Σ p log p)** calculated from source measured visit proportions among *positive* animal visitor names at each plant and original survey season, not a biologically validated effective pollen delivery rate. These same-site-context numbers do not distinguish whether anoles directly caused the difference; abundance, plant flora, camera effort, invasion history and actual predation could all differ.

**Direct common-plant/season sensitivity:** matching exact documented *plant name × calendar month* within Anijima's two natural-forest source contexts yields **10 positive source pairs**, with **5/10** showing a higher honeybee visit share in the anole-present context. **Only 3/10** have a smaller observed partner-diversity eH in the anole-present context, and the **median eH difference is +0.421 (present minus absent)**. The apparent negative crude context difference (mean eH 2.232 vs 2.719) therefore does **not transfer** to exact comparable plant-month pairs. Sampling is still very unequal (**206 vs 558 source visits within common plants**), and n=10 pairs are nested within **one** island/two observation sites; no independent-anole-site or light-mediated causal inference is warranted.

The result is a useful actual macro-ecological warning: **changes in which plants are sampled can generate apparent loss of alternative visitor channels even when the paired same-plant response differs**. This is especially pertinent to the 562-island photo-flora where the set of observed plants changes with geographic light gradients.

### Experimental ALAN and island replacement belong to separate, non-overlapping source tiers

- The **Tibet/Lu 2024 original ALAN manipulation**, source authenticated in PR #65 and common-guild measured Swiss cover comparison in #69, shows nocturnal moth/Lepidoptera visits ↓ and nocturnal fly visits ↑, but its plant-specific total-visit vulnerability class fails leave-one-block-out transfer (2/4 eligible heldout blocks replicate).
- The **2016 Switzerland/Giavi original lit/dark** shows day guild visit counts and exposure-adjusted direction strongly depend on individual field site and observed floral cover.
- The **Netherlands/Boom original 4-spectral-light experiment** has same moth as pollen mutualist and later seed predator, seed production from formula not significant net output difference.
- The **Ogasawara 2021 four-island original dataset** measures a restructured daytime community that is different in honeybee source-visit concentration and local predator records, but **no ALAN light/dark manipulation, ground lux, nocturnal census, actual moth pollination shift, stigma pollen, fruit/seed or predator feeding capture data**.
- The **2023 global 452 lizard flower-use review** distinguishes source recorded flower use versus actual confirmed pollination, but lacks accessible raw 452 interaction-level supplementary evidence for a new direct pooled island-mainland functional-efficiency estimate.

**Next defensible joined research target**: recover spatially matching local 2021 or historical plant population daylight/nighttime field cameras and effective pollen deposition, compare existing island channel state (honeybee-dominated versus diverse native/guild documented) under actual measured annual light intensity **only after the same physical island/site/date and direct effective pollen response are confirmed**. The 562-island photo/historical satellite frame is a different geographic and biological source. A source photo/light presence is not a manipulated ALAN effect.

## Reproducibility and source sensitivity

1. scripts/analysis/alan_2023_global452_lizard_flower_efficacy_original_source_probe_20261114.py — actual Wiley original source GET 403 and public separate Figshare phylogeny metadata; result artifact and CI #43575594251.
2. scripts/analysis/alan_Quitian2026_original_Ogasawara_4island_network_source_manifest_20261114.py — author original Zenodo MD5 + raw XLSX schema, CI #43575803117.
3. scripts/analysis/alan_Quitian2026_Ogasawara4island_visitation_functional_portfolio_20261114.py — raw original 2,745→2,063 positive-row source reconciliation; author visit denominator 4,174, honeybee labels, seven source contexts, 12 physical-island/month combinations, florivory versus valid visits distinguished; within-Anijima matched original plant×month sensitivity, CI #43576471698.
4. Unit tests and original GitHub Actions workflows on this Draft branch. Only aggregate source outputs, no original source plant×island×local pollinator names/coordinates copied. All earlier source results and main unchanged.
