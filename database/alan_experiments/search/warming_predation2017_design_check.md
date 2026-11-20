# Independent design check: Miller2017

## Primary-source design (reviewed 2026-11-20)

[Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5647293/) specify48 field cages in eight rows of six, with neighbouring row pairs forming four randomized blocks. Twelve combinations cross lighting, warming and predator treatments. Each illuminated cage had its own two LED lamps: the cage, not each lamp or aphid, is the lighting unit. Treatment predators were two adults of either Coccinella septempunctata or Coleomegilla maculata; all cages began with20 aphids. Aphids were collected after9days, combining field-aspirated and laboratory-counted individuals. This endpoint includes reproduction and losses, not direct consumption. Predator recovery was incomplete. The separate one-hour indoor feeding trial is a different experiment. [Original workbook](https://doi.org/10.5061/dryad.h7k26).

## Independent workbook inspection

Inspected `../gbif_annual/alan_experiment_sources/warming_predation2017/Miller data ProcB.xlsx`, sheet `Field data`, using pandas and openpyxl independently of the extraction script. The48 rows have unique `(Row,Cage)` keys; each physical row contains six cages. No workbook comments or hidden sheets provide a treatment-code dictionary.

Reconstructing adjacent pairs as `(1,2),(3,4),(5,6),(7,8)` gives four blocks. Blocks1–3 contain each factorial cell once. Block4 contains the following conflict:

| Row | Cage | Pred. | Warm | light | aphids |
|---:|---:|---|---|---|---:|
|7|4|N|Y|Y|395|
|8|4|N|Y|Y|276|

Thus46 factorial cells contain one cage, one cell contains two, and the expected block4 `CM/Y/Y` cell is empty. Across blocks, `CM/Y/Y` has3 cages and `N/Y/Y` has5; all other treatment combinations have4. The conflict is between different cages, not duplicate measurements. The article's balanced-design statement cannot identify which, if either, row is mislabelled. Neither count magnitude nor agreement with a published effect justifies reassignment.

`C7`→Coccinella septempunctata, `CM`→Coleomegilla maculata and `N`→no predator are inferred from the reported treatment identities and code names; these mappings are not an explicit workbook codebook. Preserve original codes alongside interpretations. Likewise `(Row+1)//2` is a block reconstruction consistent with the layout and first three complete blocks, rather than an original `Block` column.

## Extraction consequence

Retain all48 observations and the unresolved allocation flag. Do not force four cages per cell, replace the missing cell, drop one of the two N/Y/Y cages as a duplicate, or claim complete four-block pairing for every contrast. All C7 light comparisons have four complete block pairs, but contrasts involving warmed CM or no-predator groups require separate eligibility handling. Any complete-block subset must be labelled as a subset; preserve the full raw-data discrepancy. No numerical correction or inferential model was applied in this audit.
