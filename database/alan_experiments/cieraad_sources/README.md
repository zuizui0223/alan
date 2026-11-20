# Cieraad2023 original-data extraction

See [primary review](primary_review.md) for experimental design and sources. The original [Dataverse dataset](https://doi.org/10.34894/VD98DG), version 1.0, contains one workbook (file 320748, 60,989 bytes). Download used `https://dataverse.nl/api/access/datafile/320748`. The repository-reported SHA-1 is checked before extraction; `audit.json` records SHA-256 as well. Original workbook and metadata are cached externally, not copied here.

## Raw-data-derived results

The workbook contains 188 blueberry and 129 oak measured plants, unique IDs, 16 posts and eight illumination rows. Treatment and forest type are constant within each row. Each forest has exactly one row per treatment. No missing values occur in these retained records; this does not establish absence of losses before the analysis table was formed.

We averaged plant leaf-area-loss fractions within posts, then weighted the two posts equally within each row. The six contrasts use equal weights for the two forest blocks. Ratios of these pooled means are approximately 2.02, 4.93 and 2.64 for blueberry (green, red, white), and 0.95, 1.46 and 1.32 for oak. All approximate t1 lnRR intervals include zero. They describe retained plants at one location and are not precise population-level effects.

Forest-specific red/dark ratios are 7.77 versus 0.91 for blueberry and 2.34 versus 0.80 for oak (conifer versus deciduous). With one row per treatment in each forest, these contrasts have no independently estimated within-forest precision. They are descriptive heterogeneity, not replicated forest-type tests or moisture effects.

`effect_covariance.csv` retains shared dark rows, species and forest blocks. Its rank is one: six contrasts do not supply six independent experiments. Main catalogue inclusion remains held for programme overlap and final adjudication. Plant-level published regression precision is not substituted for row-level precision. Growth, gall counts and leaf traits have not been extracted as meta-analysis effects here.

Reproduce with `python database/alan_experiments/cieraad_sources/extract.py`. All source cell counts and row means are retained for checking.
