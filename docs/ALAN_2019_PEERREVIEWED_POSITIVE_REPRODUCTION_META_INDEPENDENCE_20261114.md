# ALAN 2017–2026 reproductive response reversals: an independent peer-reviewed positive counterexample and the *actual* source-data gate

**Study (new to current ALAN source atlas):** Macgregor CJ, Pocock MJO, Fox R, Evans DM (2019), *Ecosphere* 10:e02550, DOI [10.1002/ecs2.2550](https://doi.org/10.1002/ecs2.2550). This paper was published **18 January 2019** (not a new 2026 result) and is an independent **peer-reviewed experimental study**, separate from Li et al. 2026 preprint and the Knop et al. 2017 Swiss negative plant-fitness result.

## 1. Why this changes the project hypothesis

Original Macgregor experiment uses ***Silene latifolia*** **female** plants surrounded by male pollen-donor plants, with each female group allocated to full-exclusion, day-only, night-only and open pollinator-access cages. Female plant success depends on pollen transfer: its positive response is **not rescue through autonomous selfing**. Two lamp types (high-pressure sodium vs LED) and two timing regimes (all-night vs until midnight) were compared with darkness in the original experiment. Original paper **found a higher probability of a flower being successfully pollinated in FULL-NIGHT lit treatment than in dark control**, and **no clear corresponding improvement under PART-NIGHT light**. **Seed count and total seed mass per pollinated capsule did NOT show a resolved full-night-vs-dark increase** (different endpoints). Direct individual insect movement, effective pollen contribution by insect guild and true net seed production are not identified from the reported odds of per-flower success. This belongs to the primary authors.

This genuinely **falsifies the overly strong working hypothesis**: "plants unable to self or with nocturnal visitors ALWAYS lose reproduction under ALAN." It also **does not prove** that the cause is added moth visits or beneficial light corridors. The study notes lit plots were powered by petrol generators absent from unlit controls, and possible other indirect environmental confounding. Its light treatments and pollinator access allow conditional analysis but one study cannot establish a universal positive effect.

## 2. A non-cherry-picked counterexample framework

| Primary published paper | Distinct evidence unit | Light contrast | Reproductive endpoint | Original source result | Evidence status |
|---|---|---|---|---|---|
| Knop et al. 2017 Nature doi:10.1038/nature23288 | One Swiss manipulated community; focal Cirsium | illuminated/dark | focal fruit set | **−13%** (published) | peer-reviewed; only authors' published effect represented |
| Macgregor et al. 2019 Ecosphere doi:10.1002/ecs2.2550 | One UK manipulated plant-caging study, four experimental time runs | full-night/dark vs part-night/dark, HPS and LED | per-flower successful pollination | **positive under full-night**, part-night unresolved; no resolved seed-count increase | peer-reviewed published conclusion, NOT source row reanalysis |
| Boom et al. 2020 Sci Rep doi:10.1038/s41598-020-75471-1 | One Dutch four-field S. latifolia–Hadena study | white/green/red/dark spectrum | fruit set, larval infestation, ovules and estimated seeds per flower | strong constituent changes; **net seeds ANOVA p=0.63** | author table reproduced in ALAN #67 |
| Li et al. 2026 bioRxiv doi:10.64898/2026.07.10.737511 | One UK three-hayfield/79-plant source study | 26-bulb array / 4-bulb point / dark | any seed per plant, fruits and seeds per plant | array vs dark **20/27 versus 8/26**; positive 3/3 fields, point inconsistent | **unreviewed preprint** actual author raw CSV reanalysis in ALAN #74 |
| Lu et al. 2026 Journal Plant Ecology doi:10.1093/jpe/rtag160 | One Tibet two-botanical-species experiment | ALAN/dark × four pollinator access arms | stigma conspecific pollen, seed set | *Pedicularis kansuensis* **−56% seeds**; *Euphrasia regelii* no resolved seed effect | peer-reviewed accepted article abstract; **no original rows refitted** |
| Lu et al. 2026 Ecology & Evolution doi:10.1002/ece3.74368 | Tibet 2024 five-block animal visitation network | ALAN/dark × day/night | **insect VISITS, not seeds** | nocturnal Lepidoptera **−57.4% visits**, flies **+179.6% visits** | original source 400 rows, 5 fields replayed ALAN #72 |

Six publications, only **five independent author-research-group labels**; the *two Lu 2026 papers have four authors in common*, an alpine Tibetan ecological setting and potentially same broader project. In Ecol Evol source, actual site = Gansu Gannan National Research Station, **summer 2024, five paired blocks, 8000 K LED**. The JPE 2026 accepted-manuscript abstract does **NOT disclose exact source plot IDs/date**, so neither true overlap nor full independence can be certified. **Do not treat them as two independent illumination experiments in a meta-analysis** until original field design matches are audited.

Source records also separate **fruit set**, **per-flower pollination success**, **per-plant ≥1 seed**, **absolute seed count**, **stigma pollen**, and **visitation**. These have different sampling units; simply converting published percentages to a shared slope and adding a random-effects model would not be a valid new meta-analysis.

## 3. Crucial authentic source inspection — original R code is accessible, original numerical CSV isn't in these verified public archives

The original Macgregor article links [Zenodo DOI 10.5281/zenodo.1743917](https://doi.org/10.5281/zenodo.1743917). Verified GitHub Actions REST record metadata **ID 1743917**, 1 ZIP, **16,180 bytes** with declared MD5 **26c45205c4bcef4576be592d54893615**. The actual authentic Zenodo ZIP was downloaded and its members enumerated: it contains only six \`.R\` scripts (no raw original \`.csv\` files).

The author's public repository [CallumJMacgregor/BioassayPaper](https://github.com/CallumJMacgregor/BioassayPaper) was independently checked. SHA-verified original \`FinalScripts.R\` Git blob **04136ea82d6c1bbef22f4fd60c41de9907b00550** calls:
- \`Data/SeedSetBinom.csv\`, for original \`cbind(Successes,Failures)\` binomial GLMM;
- \`Data/SeedSetYN.csv\`, for 0/1 per-flower success and Fig. 4 full-vs-part-night light predictions;
- \`Data/SeedCount.csv\`, for actual seeds/capsule and seed mass analysis.

**All three necessary original experiment CSVs are ABSENT from the independently verified author GitHub tree and Zenodo archival ZIP.** There may be other repositories or author-held data not checked. The archived script uses \`(1|fRound)\` as a random effect, adjusts pollinator access, plant usage and lamp distance in \`model2YN\` and \`model3YN\`; no new coefficients, SE, sample counts, held-out-plot test or causal effect were numerically computed here. A code archive alone is NOT a source reanalysis.

Original source verification: [Actions #43583533873 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43583533873) with seven synthetic/identity tests, original author \`.R\` blob SHA and actual Zenodo ZIP member listing; compact signed artifact **#14108984434**, ZIP digest \`sha256:3982298519cfdeab8609af3219164c2ff2cf83c3f65bf25c725b8be284d3f4a1\`. Earliest Zenodo metadata-first [Actions #43583375445 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43583375445).

## 4. Revised primary scientific hypothesis

**Not** "night-pollinator-reliant or self-incompatible plants always decline with ALAN."

Test instead: *net reproductive susceptibility* is the interaction of:
- source-verified original nocturnal share of effective pollen delivery, AND verified autonomous selfing or other genuine reproductive assurance;
- local lighting **exposure** regime (lux, spectrum, duration, spatial array, potentially lamp count);
- whether visitors are displaced **from flowers** versus spatially aggregated **near lit flowers**, and whether increased non-moth visits are truly effective;
- flower availability and diel/seasonal overlap; negative adult moth predation and positive flower-pollen transfer by nonstandard visitors;
- plant population history, island type, pre-existing true pollinator shifts and actual functional-channel redundancy.

The strongest new falsification is **a peer-reviewed positive reproductive response even in a system using female flowers that cannot self**. Thus selfing is a *potential buffering moderator*, not the sole switch in response sign. The same Silene genus/species also participates in nursery pollination/larval seed predation under spectral lighting (Boom et al. 2020), so a fixed plant species name is not enough to predict response without local visitor services.

## 5. What to do next — inspectable stopping gates

1. Retrieve **Macgregor 2019 original THREE CSVs** from the authors or other legitimately shared data archive, then rerun actual \`model2YN/model3YN\` per-flower GLMM and add whole-experiment-run / plot sensitivity. Until then use only the **published direction** and no self-fitted variance.
2. Retrieve **Lu JPE 2026 primary field sites, years, block identifiers and original reproductive row-level source** (currently original published abstract/DOI only), match against Lu Ecol Evol 2024 five fields before assuming independent replication. No fabricated sample overlap verdict.
3. Conduct a true outcome-specific meta-analysis of **fruit/pollination success (binomial per flower)** only when both treatment×control counts/SE and independent experiment units are established. Analyze seed numbers **separately**, not pooled with visits or fruit success. Group all outcomes by actual original experiment and field to avoid pseudo-replication.
4. Evaluate *external moderators* with sufficient independent study groups and levels: spectrum, duration, array geometry, selfing, true night dependence. At present 6 papers / 5 group IDs and outcome heterogeneity are insufficient for a statistically supported global 4-way interaction.

## Reproducibility

New code and tests:
- \`scripts/analysis/alan_Macgregor2019_original_zenodo_and_study_independence_20261114.py\`
- \`scripts/analysis/alan_Macgregor2019_original_script_data_dependency_20261114.py\`
- \`tests/test_alan_Macgregor2019_original_study_independence_20261114.py\`
- \`tests/test_alan_Macgregor2019_original_script_data_dependency_20261114.py\`
- \`.github/workflows/alan-Macgregor2019-metaanalysis-study-independence-20261114.yml\`

Only source archive member names, code identity/expected CSV filenames and publication-level effect tiers exported. No author raw flower records, author code, dataset binaries or species×original rare-site observations redistributed in this repository. Keep Draft stacked on previous ALAN #74 and main unchanged.
