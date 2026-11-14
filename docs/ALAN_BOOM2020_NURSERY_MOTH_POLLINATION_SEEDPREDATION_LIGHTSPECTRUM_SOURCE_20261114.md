# ALAN — the same moth is a pollinator and seed predator under spectral nightlight

**Independent mechanistic field study:** Boom et al. (2020), Scientific Reports, [DOI 10.1038/s41598-020-75471-1](https://doi.org/10.1038/s41598-020-75471-1), [PMCID PMC7591485](https://pmc.ncbi.nlm.nih.gov/articles/PMC7591485/). One original experiment across **four Netherlands field sites**, not four globally independent published studies. The biology and original results were already published by Boom et al.; our contribution is an independent source-locked audit and a clearly bounded net-service decomposition.

**Real primary source replay SUCCESS:** [GitHub Actions #43570186936](https://github.com/zuizui0223/alan/actions/runs/43570186936). The script downloaded actual original Europe PMC JATS XML, confirmed journal DOI and all four rows/flower, fruit and ovule nested denominators in printed Table 1, and recreated the original estimated reproductive-output formula. The original public JATS receipt SHA256 is 182dafe1e99fef5453ecf0d1d1401361fe949fd8d8617ae040588fd19795154f. Experimental light was ground-level approximately **7.6 ± 1.2 lx**, not satellite PCNL. Treatments: darkness, white, green and red lighting.

## The core functional conflict

The moth *Hadena bicruris* pollinates *Silene latifolia* flowers, and its larvae consume seeds inside fruits following female oviposition. Thus the same taxon is a **mutualistic pollen vector and antagonistic offspring seed consumer**. This is distinct from gecko or bat predation on adult moths and potential direct pollen-service substitution, neither of which was measured in this experiment.

| Published original treatment | Flowers producing fruits | Fruits infested | Fertilized ovules among sampled noninfested fruits | Published estimated seeds per flower |
|---|---:|---:|---:|---:|
| Dark | 0.73 | 0.28 | 0.56 | 117 |
| White | 0.78 | 0.44 | 0.49 | 86 |
| Green | 0.75 | 0.53 | 0.43 | 61 |
| Red | 0.71 | 0.33 | 0.47 | 90 |

Original inferential conclusions: fraction fruit-setting **p=.17, unresolved**; fruit infestation varies by light spectrum **p<.001**, ovule fertilization **p<.01**, but original estimated seeds-per-flower ANOVA **p=.63, unresolved**. No meaningful significant net seed loss is inferred. Original seed outputs are model-derived from treatment means, not individual mature seeds counted in an independently replicated global trial. Flowers, fruits and ovules are **nested within plant, transect and site**, and the original paper accounts for this.

## Descriptive formula audit

The authors' original algebraic seed-output definition: seeds per flower = p(flower sets fruit) × [1 − p(fruit infested)] × p(ovule fertilized in selected noninfested fruit) × mean ovules per fruit (404 ± 7.78).

Exact arithmetic on original *rounded* treatment means:

| Term | Dark | Green | Green/dark factor |
|---|---:|---:|---:|
| Fruit set | .73 | .75 | 1.027 |
| Escaping Hadena infestation | .72 | .47 | 0.653 |
| Ovule fertilization among uninfested fruits | .56 | .43 | 0.768 |
| Estimated seeds/flower from rounded components | 118.91 | 61.24 | **0.515** |

The product is **−48.50%** green/dark, a reproduction of the authors' *arithmetic*, not a newly significant field result. Published displayed means were **117 vs 61** because unrounded quantities differ slightly. On the log product scale, factor contributions are fruit set +0.027, seed-predation escape −0.427, sampled uninfested-ovule fertilization −0.264. These are identities, **not a regression mediation decomposition** or causal attribution: ovules were sampled conditional on noninfestation, which can cause selection bias.

## Implications for syndrome vulnerability

The primary ALAN estimand must distinguish **actual pollen deposition, fruit set, net viable seed output and seed loss due to the visitor's offspring**, in addition to flower visitation. A nocturnal specialist may maintain pollination yet shift toward antagonistic seed consumption under particular light spectra. This makes a simple moth-visit-loss proxy an incomplete syndrome-vulnerability criterion.

Contrast with independent original systems, keeping endpoints separate:
- Tibetan meadow Lu 2024: illuminated nocturnal Lepidoptera visits **838→357**, Diptera **54→151**; original 11 plant-group pooled results are **not robust** to held-out blocks.
- Swiss meadow Giavi 2016: daytime Diptera count **569→550**, despite positive within-pair direction in **5 of 6** sites, and wide whole-site bootstrap interval; plant/guild daytime effects not transferable as a single universal sign.
- Netherlands Boom: pollen-vector and seed-consuming **same moth**, original fruit/infestation/ovule components across 4 light spectra with nested samples. Neither ALAN experiment alone resolves island-vs-mainland shifts, predator replacement service or true local phenological mismatch.

**Not established:** a global ALAN-by-nursery-pollination coefficient, effective replacement of moths by geckos/bats/centipedes, adult moth predation pressure, a nightlight-caused phylogenetic shift on islands, a net significant seed reduction in this Boom trial, or a mismatch between flower opening and the actual local pollinator schedule.

## Reproducibility

Source original parser: scripts/analysis/alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114.py. Tests: tests/test_alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114.py. Workflow: .github/workflows/alan-original-Boom2020-nursery-dual-role-spectral-20261114.yml. [CI #43570186936](https://github.com/zuizui0223/alan/actions/runs/43570186936) SUCCESS; only public source treatment-level aggregates in the Actions artifact. Draft stacked on #66; older PRs and main unchanged.
