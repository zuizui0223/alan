# Bat forest-edge refuge: raw-data support audit

2026-11-21. [Primary article](https://doi.org/10.1186/s40462-020-00238-2); [author data and code](https://doi.org/10.5281/zenodo.4036279), CC BY 4.0. Primary Methods/Results were rechecked from the cached article XML; the author analysis code was read, not executed. This is the existing Dutch Licht op Natuur programme, not an independent programme added for each paper or bat group.

## Original-data reconciliation

The public CSV contains 35,067 positions. Reproducing the author's selection of trajectory identifiers containing more than one position retains **28,646 positions in 1,855 trajectories**. Every retained position has recorded spatial imprecision at most 1 m. The group totals exactly reproduce the article: Eptesicus/Nyctalus 1,788, Myotis/Plecotus 692, Pipistrellus 26,166. Trajectory IDs each map to only one site, night, treatment and group. A trajectory is not an independently identified animal; repeat passages remain possible.

| Treatment | Sampled nights | Unique site–treatment cells |
|---|---:|---:|
| Control | 4 | 4 |
| Red | 4 | 3 |
| White | 4 | 3 |

There are **4 sites, 12 nights and 10 site–treatment cells**, not four independent sampled sites per colour. Site 67 contributes red observations on 10 and 21 July 2018, and white on 12 and 22 July. Sites 1 and 67 have all three conditions; site 3 has control/white, and site 45 control/red. Treatments were recorded on different nights, not simultaneously. The original Methods assign one row per colour within site, so the repeated site–colour observations are not new treatment allocations. Raw recording IDs differ and are retained in `sampled_nights.csv`.

`night_group_descriptive.csv` retains observed night–group cells with counts, distance and forest-position fraction. Missing night–group combinations are not filled with zero. These fractions are descriptive summaries of detected positions, not distance-standardized effects, independent animal proportions, total bat abundance, feeding success, or predation rate.

## Why no pooled predation effect is added

The measured response is flight position relative to the forest edge. The original models describe a treatment-by-distance response: forest use increases near lights in group-dependent ways. Exponentiating a light coefficient alone would evaluate a distance-specific contrast at the reference distance, not a general lighting effect. Coefficient covariance and common distance support are needed for any defined prediction contrast.

The published model table and saved code also require reconciliation before numerical reuse:

- For Eptesicus/Nyctalus, all-position white×distance SE is printed as **0.435** in Table 1 but **0.04347** in the saved code output (estimate approximately −1.639 in both). Neither value is silently corrected.
- For Myotis/Plecotus, the table's all-position distance coefficient **−1.864** matches the code output with white as reference. The code with control as reference gives **−0.3277**. Thus the displayed base slope cannot be uncritically combined with control-reference treatment coefficients.
- The code uses precision weights in the Eptesicus/Nyctalus full model; the Pipistrellus and Myotis/Plecotus full models shown are unweighted. These are not three identically weighted analyses.

The author models include recording/night ID and trajectory random intercepts, but do not thereby establish independence of individual bats or separate every site, treatment and date influence. No replacement model or variance based on 28,646 independent positions is fitted here. No feeding-buzz, prey identification, attack, or pollination endpoint is supplied by this response.

This recovery adds zero effects to the 108-row pooled-effect catalogue and does not inflate the independently scoped reported-statistics registry. It provides an auditable alternative predator response—spatial avoidance/refuge use—rather than evidence that lighting universally raises predation pressure.

## Reproduction and provenance

Download the two files listed in `audit.json` and save the Zenodo API record as `record.json` in a source directory. Run `python database/alan_experiments/bat_refuge2021_sources/extract_support.py SOURCE_DIRECTORY`. It checks API byte sizes and MD5, records SHA-256, reconstructs the source selection and writes the three descriptive tables. Raw files are cached outside the repository; the 25 MB position table is not duplicated here.
