# Independent UK NPMS test of ALAN × plant moisture affinity (2026-11-17)

**Status: original-source empirical fit executed and signed, observational, post-hoc extension.** The original alan/main finding (283 nightlight-associated recording decline candidates, 158 moisture-supported cases, environment-adjusted EIV difference +0.53384) remains **unchanged**. This separate validation does not use the 283-case selection to choose its plant species or outcomes.

### Original independent sample

- Source: [UK National Plant Monitoring Scheme 2015–2025](https://doi.org/10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8); original repeated positive-cover field-data source from [ALAN run #42815428180](https://github.com/zuizui0223/alan/actions/runs/42815428180) and original annual 2014–2024 historical Li–Cao [PCNL v3](https://doi.org/10.5281/zenodo.17013414). The satellite observation precedes the botanical survey year. This is not ground-level illumination.
- The prior source was independently verified at exact raw positive-cover SHA256 **4e8bb59f3704d0c6140194fdbf0a3d1f50e11235da243b894833b21824d3a165**; frozen original year-lagged PCNL source SHA **5ecae6a9b62de636c8e68fd694c69698eec154e5a8a0ba641beb5a7e6061e427**.
- Exact independent [Tichý et al. (2023) European EIV](https://doi.org/10.1111/jvs.13168) with unchanged public CSV SHA256 **ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba**. Only exact species binomials accepted, no guessed aliases or unmeasured local moisture value assignment.
- Before EIV joining: **61,555** positive-cover plant×plot×method×season×year values, **686** botanical species, **1,138** independently identified original plots. Eligibility is at least two positive-cover record years for the *same* plot/species/method/season and ≥0.02 log1p previous-year PCNL span.
- After exact MOIST species matching: **30,056** positive-cover rows / **438** species (**143** MOIST>7, **295** MOIST≤7) / **1,054** plots.
- Strict shared-local-environment panel: **2,178** truly jointly observed wet/dry plot×year×method×season cells in **389** unique plots, across **11 calendar years**.

### Independent source result (not ALAN causality)

Outcome is the original DOMIN-coded *recorded positive cover* midpoint, converted to percentage points. It is **not** total plant abundance, local disappearance, failed recording of a species, flowering or realized reproductive success. Original source survey ID and broad season are held constant inside panel. Main estimates use physical plot-equal weights, plot-cluster nominal uncertainty and Tikhonov-free weighted sparse FWL fixed-effect residualization.

| Contrast | Design and blocking | Coefficient on log1p prior PCNL × MOIST | Nominal plot-cluster 95% interval | p |
|---|---|---:|---:|---:|
| Primary source-only continuous wetness interaction | plot×species×survey×season FE + species×year FE + source survey day adjustment | **+0.226659** | **[−0.105978,+0.559295]** | .1815 |
| Stricter same-plot-year wet/dry environment | above + plot×year×survey×season FE (absorbs main light) | **+0.293302** | **[−0.087962,+0.674567]** | .1312 |
| Wet >7 binary source moderator | primary FE | +0.344692 | [−0.859649,+1.549033] | .5745 |
| Wet >7 binary strict within plot×year | strict FE | +0.500894 | [−0.687595,+1.689382] | .4078 |

Units of the *continuous interaction* are **recorded cover percentage points per one EIV MOIST unit × one unit log1p annual prior PCNL**. Binary interaction units compare wet (>7) and drier (≤7) plants in prior-PCNL slope. None of the four contrast intervals excludes zero. Their point estimates are **positive**, opposite the simplistic hypothesis that higher wetness entails a more negative ALAN slope of *positive cover*. This is not evidence of a positive causal ALAN effect either.

In the primary continuous model, the residual *main light* slope at MOIST=7 is +0.7501 [0.0105,+1.4897], nominal p=.0468. This one nominal coefficient is part of a previously selected trait-covered subcohort and **must not be promoted as a plant-benefit claim**; the previous full 686-species NPMS source cover model had +0.370 [−0.095,+0.835] with plot-equal weights, while record-equal and panel-equal estimates were near zero. The strict model absorbs common site-year light effects.

### Scope and interpretation

This new result constrains one candidate mechanism: a broad **ALAN × inherent moist-habitat affinity response in recorded plant cover** was not confirmed. It does not falsify the main alan/main finding that the **selected long-term recording-decline species** have higher wetness affinity. Cover *conditional on positive recording* cannot measure extinction or unreported plants, and different sampling systems and observation years can give different ecological questions.

Other mechanisms remain plausible without assigned contributions: habitat hydroperiod/drainage and mowing changes correlated with anthropogenic development; ALAN×nutrient competition/plant phenology; visitor behavior and reproductive timing; water-derived insect supply and predator/herbivore response. The wetness proxy is not local water table and the satellite radiance is not leaf lux. Observer season/protocol and shared spatial climate/urban trends remain confounding threats. The study is a geographically UK-weighted survey, not an independent global replicate.

**Next preregistered stress:** original NPMS source light-variation cutoffs **≥0.1** and **≥0.3** log1p prior PCNL to check whether small annual satellite-light changes cause misleading signs. These should be reported regardless of direction, with no coefficient if source common support fails.

### Reproducibility and provenance

- [Original numerical fits and conservative site/trait-coverage gates](../results/all_species_forest/mechanism_candidates/npms_moisture_original_repeated_positive_cover_20261117.json).
- [Source rehydration, original hash validation and fixed-effect model](../scripts/analysis/alan_npms_independent_moisture_priorlight_20261117.py).
- [Synthetic exact-coefficient tests](../tests/test_alan_npms_independent_moisture_priorlight_20261117.py).
- [Original Actions #44074149632](https://github.com/zuizui0223/alan/actions/runs/44074149632), both workflow jobs success, including all three numerical synthetic checks; original NPMS source SHA and original European EIV source SHA validated. No raw 30,056 location×plant individual records committed.
