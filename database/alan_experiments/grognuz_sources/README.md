# Floral resource timing: original-data version audit

Article: https://doi.org/10.1016/j.biocon.2025.111650 . Retrieved data: https://doi.org/10.17632/xjmn5jhxyn.1 . Article-cited version: 10.17632/xjmn5jhxyn.2 .

## Primary-study scope

Grognuz, Gisler and Knop studied Swiss wildflower strips in 2022–2023, using six lit–dark pairs in 2022 and eight in 2023. These are site-year counts, not proven distinct sites across years. Germination used 2023 sites; the article reports 22 sown species and six selected for analysis. Phenology and petal movements measure resource timing, not observed plant–pollinator temporal mismatch. The authors report no clear germination-probability effect and species-dependent timing effects. Signed timing differences require day/minute effect scales, not a uniform lnRR of dates. Sources: Methods2.1–2.5 and Results3.1–3.3.

## Acquisition audit

The version2 landing page returned HTTP200 but its embedded dataset state returned404/not found. Version1 was accessible. Three small original files were downloaded and matched their provider SHA256 values. The 55MB flowering table was inventoried but not downloaded in this step. The author-institution PDF was obtained and hashed. See `source_manifest.csv`; no claim that all original data were retrieved.

The version1 germination workbook contains 6336 rows, 11 species codes, 16 sites and 36 observation days (DOY132–167 in2023, May12–June16). It has unique site/species/day keys. The cumulative germination field does not decrease. Repeated cumulative values must not be summed as daily germination events. Article methods describe sowing May10 and monitoring through June20; do not assume the available window equals the final reported window.

At the latest available day, nine species codes meet the stated minimum of at least five germinations in each light arm and germination at at least three sites per arm. The paper describes six selected species. This discrepancy needs version/filter reconciliation; it is not evidence that either dataset is erroneous. `germination_v1_support.csv` and `germination_v1_last_available.csv` preserve the check without inventing names from abbreviations or treating absent species as zero.

All16 germination site labels match the original2023 site list in the separately audited slug/herbivory dataset. This supports a shared-programme flag, not proof of identical individual plants, outcomes or treatment exposure. Site-year treatment and coordinates need further verification before covariance modelling. See `slug_programme_site_crosswalk.csv`.

No effect sizes are added to the main catalogue. The version discrepancy and potential shared programme remain explicit holds. Do not present the paper's model estimates as independently reproduced, or infer pollinator mismatch from flowering changes alone.
