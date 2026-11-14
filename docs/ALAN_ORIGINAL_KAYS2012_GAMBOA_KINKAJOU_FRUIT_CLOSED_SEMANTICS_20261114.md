# ALAN — real Kays 2012 balsa study links original same-flower visitors, fruit production, and an incomplete control

**Original author/source replay SUCCESS (2026-11-14):** [GitHub Actions #43562922775](https://github.com/zuizui0223/alan/actions/runs/43562922775), [aggregate artifact #14101138236](https://github.com/zuizui0223/alan/actions/runs/43562922775/artifacts/14101138236). This is stacked on [Draft PR #62](https://github.com/zuizui0223/alan/pull/62); main and all previous branches remain unchanged.

## Same published experiment / year / area — independently verified

The original **Kays et al. (2012), Mesoamericana 16(3), pp. 56–70**, [public published original article (PDF)](https://www.revistas.up.ac.pa/index.php/mesoamericana/article/download/907/777), reports animal visiting behavior, phenology and four flower-fruiting manipulations on balsa (*Ochroma pyramidale*) trees in central Panama. The 2009–2010 flowering season includes the original Kendall/Nicholson meta-analysis **2010 January–February** fruit-set treatments. The authors' original 1652-row source contains `Kays_2012` with exactly four original arm means and sample sizes. The original meta-analysis coordinate of the Kays fruit-set experiment is **35.58 metres** from the paper's Gamboa experimental tree/tower, and approximately **3.285 km** from its Pipeline tree and **15.98 km** from its BCI tree. Only distances are exported to the aggregate science result; no original plant/site-year coordinate rows are replicated in this repository.

Kays source is therefore a **second independent original published study**, after Wenzell 2024, whose real visitor behavior and fruit reproductive treatments can be anchored to the same original field area and flowering season. Unlike Wenzell's population-level observed visitors, Kays' Table 3 describes animals **observed on particular experimental flowers**. Not all visitors could be continuously observed at night, so visitor identity and effective pollen donation are not experimentally isolated.

## Author-printed reproductive experiment rehydrated exactly

| Original treatment | Actual fruited / source flowers | Fraction | Protocol meaning |
|---|---:|---:|---|
| Day-only open | 1 / 13 | 0.076923 | flowers exposed to diurnal animals, netted after sunset |
| Night-only open | 3 / 6 | 0.500000 | flowers unnetted after sunset, other periods excluded |
| Open both | 4 / 6 | 0.666667 | unrestricted |
| Wire-mesh "small insects only" | 0 / 13 | 0 | **bees can pass** through ~0.9 cm mesh; NOT complete exclusion of animal visitors |

All original author source treatment means and `sample_size` values match published primary Table 3 exactly, in the original source with Git blob **39a82d966afcdde670603add1c881ace8791ffe9**. The original 301 study×site×year source is independently verified to **SHA256 21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af**. Same four-arm source row and original 153/79/85 full cohort are reconstituted, without changing original labels.

**Actual behavioral evidence at original experimental flower grain:** the primary article reports **kinkajous feeding on all three night-only exposed flowers that subsequently fruited**. Other mammalian flower feeders can also visit; observations did not extend through the full night, and this association cannot prove kinkajou-specific pollen deposition or exclude bats/other animals as co-vectors. Kays describes pollen carried visibly on arboreal mammal faces and animal predation/competition around flowers, including mantises and snakes near nectar-seeking animals. This is ecological multi-channel flower-community evidence, **not** an ALAN treatment result.

## Important misclassified functional control

The original publication explicitly explains that **bees passed through its wire-mesh cage**, while vertebrates and larger moths could not. The global author source nevertheless labels its numerical 0/13 arm `complete exclusion`. Its **zero fruit set is real**, but the manipulation is not a full no-visitors bagging control and cannot be read as demonstrated autonomous-selfing absence; it is a restricted visitor-access outcome with small insects still possible.

The treatment labels `day pollination` and `night pollination` are genuine temporal access manipulations, but **never isolate one moth, one mammal or a fully monitored animal community**. Protocol audit must precede cross-paper factorial interpretation.

## Refit original global 153 four-arm cohort with the one misclassified control excluded

Independent study medians (not plant-year rows treated as independent) were recomputed from the original published treatment source:

| Original experimental summary | Full 153 source quartets / 79 papers | Excluding Kays, 152 quartets / 78 papers |
|---|---:|---:|
| Four-arm I: (O+C−D−N)/(O+C+D+N), study median | −0.05747065 | −0.06041627 |
| Original independent paper medians with I negative | 49 | 49 |
| I two-sided nominal study sign p | .04216544 | .03078650 |
| I Holm-adjusted across two newly explored endpoints | **.08433089** | **.06157299** |
| J unrestricted versus best one-window study median | +.01449544 | +.01415110 |
| J two-sided nominal p | .11001560 | .13538511 |
| K open versus "closed" study median | +.83793738 | +.80720399 |

**Scientific interpretation**: the earlier global recurrent subadditivity trend remains negative after the scientifically motivated Kays omission, but multiplicity-adjusted p still exceeds .05. **The omission does not establish a global, causal daytime/nocturnal visitor substitution mechanism.** It removes one semantically problematic experimental control and is one highly overlapping source sensitivity, not another independent global replication. The basic open-vs-closed association is already in original Kendall/Nicholson source literature.

## Independent biological guild evidence count

- **Wenzell et al. (2024)**: one independent paper with six matched 2019 plant population-years (five strict <=1km, one provisional ~1.65km); hawkmoths/bee/hummingbird activity, direct fruit-set outcomes. Previously published ecology, original-data bridge in PR #62.
- **Kays et al. (2012)**: second original independent paper with Gamboa experimental flowers, Jan–Feb 2010 in original 2009–2010 flowering season, **~36m Gamboa coordinate match** and observed kinkajou access to each of three fruiting night-only flowers. Primary source now verified here.
- **Marques et al. (2015)**: third independent article with actual bat/hummingbird direct fruit/seed experiments in same study area, but visitor-event year not independently certified; remains *paper-area evidence*, not a strict complete experiment-year bridge.
- **Ganesh & Devy (2000)**: another original predator-pollinator paper with 1996–1998 study and original source year 1998, not yet an independently verified individual source-year visitor concordance. Different biological effect measures and possible physical damage, no source matched taxon-specific fecundity.

Treat this as progress in **validated study-level source coverage**, not a sample size of nine independent studies or 25 macro-confirmatory papers. Papers with some observed fauna do not automatically offer same-site/year identity, and direct visits do not equal effective stigma pollen.

## Source-code documentation

Analysis: `scripts/analysis/alan_original_Kays2012_sameflower_vertebrate_protocol_20261114.py`; tests: `tests/test_alan_original_Kays2012_sameflower_vertebrate_protocol_20261114.py`; workflow: `.github/workflows/alan-Kays2012-original-flower-guild-protocol-20261114.yml`. Original published experiment source and original SHA+blob gates passed in [#43562922775](https://github.com/zuizui0223/alan/actions/runs/43562922775). Aggregate-only JSON [artifact #14101138236](https://github.com/zuizui0223/alan/actions/runs/43562922775/artifacts/14101138236). The original source *already includes* Kays; no artificial eighth meta source study or rewritten historical ALAN slope.
