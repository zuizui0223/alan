# ALAN B | Full original 340-reproduction panel: stage-dependence universality FAILS (2026-11-15)

## What was newly executed

SHA256-authenticated Kendall & Nicholson (2025, *Ecology Letters*, DOI 10.1111/ele.70036) author data from original pinned GitHub commit `1945b7c114ad87000608ce6018d0280001c04ef1`; **1,652 source treatment rows → 340 reproductive day/open panels in 113 original paper IDs**. Unlike earlier historical **42-row** ALAN extraction, the full source match requires identical paper, plant/cultivar, country, locality, coordinates, exact site period, start/end years/months and additional experimental treatment. It refuses ambiguous multiple stage cells or impossible fruit probabilities >1.

The signed GitHub Actions [#43712772688 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43712772688) ran the original author source re-download and **7/7 tests**. [Source script](../scripts/analysis/alan_B_global340_paired_stages_20261115.py), [tests](../tests/test_alan_B_global340_paired_stages_20261115.py), [permanent source-derived summary](../results/alan_B_global340_source_stage_summary_20261115.json). An ephemeral full source-pair artifact also exists on this run; no author's original rows redistributed.

## True new result: prior 6-paper stage tendency does NOT generalize

Matched original fruit probability to reported seed endpoints: **54 exact pairs** (coarser key 26 pairs is **audit-only**, not merged because oversimplifying keys collapses distinct sites/times and creates ambiguities). No missing arms set to zero.

| Seed endpoint matched to fruit formation | Exact pairs | Source papers | Paper-equal median of (V_seed−V_fruit) | Direction / uncertainty |
|---|---:|---:|---:|---|
| Seeds per **formed fruit** | **49** | **18** | **+0.01103** | 9 positive, 8 negative, 1 zero; paper-resampling *descriptive* 95% interval **[-0.393, 0.101]** |
| Seeds per original ovule among fruit | **1** | **1** | **+0.358** | Single Nakajima 2018 paper, ungeneralizable |
| Seeds per plant | **4** | **3** | **+0.070** | 3/3 source papers positive, posthoc sign two-sided probability 0.25; too sparse |
| Seeds per initial flower / seed-formation probability | 0 | 0 | unavailable | Cannot infer total original-flower seed success |

For the broadest *comparable conditional seed-per-fruit endpoint* (49 source pairs, 18 papers), the paper-level point-median difference is ~zero. The paper-level sign diagnostic gives two-sided probability **1** (not a registered test, nor an original-unit variance-corrected p). **There is no evidence in this panel that nighttime-access benefit is systematically stronger at the conditional seed-than-fruit stage.** The contrast is descriptive, not an effect of nocturnal artificial lighting.

### Why the historical 10/10 positive contrasts are not a general finding

The previously verified historical 42-row extraction retained **10 strict source pairs / 6 papers**, all showing V_seed > V_fruit. That table draws from **different source papers and selection decisions**, not a nested subset of the full Kendall author 113-paper reproduction panel. Its six paper IDs were `Stone_2018`, `Morse_1983`, `Robertson_2021`, `Jabbari_2013`, `Vanderplank_2020`, `Noone_2022`; the new full-source matched endpoint set shares just **2/6** source paper IDs: `Morse_1983`, `Vanderplank_2020`. Hence the old 10/10 result remains factually valid *for those selected source rows* but MUST NOT be called a replicated global gradient or pooled with the 54 pairs as independent data. Differences can also arise from seed endpoint definition and digitization/subset rules.

## All 18 matched source-paper seed-per-formed-fruit stage deltas

V = 1 − (DayOnly/Open). Negative deltas remain **negative**, not removed. Multiple stratum effects within one original paper were reduced to its median before counting directions or paper-resampling.

| Source original paper ID | Median V_seed/fruit − V_fruit |
|---|---:|
| Dar_2006 | -1.497 |
| Borbon-Palomares_2018 | -0.813 |
| Cruz-Neto_2015 | -0.557 |
| de Avila_2011 | -0.433 |
| Hove_2007 | -0.424 |
| Alonso-Pedano_2012 | -0.392 |
| Larrea-Alcázar_2011 | -0.034 |
| Luo_2011 | -0.024 |
| Casas_1999 | +0.000 |
| Aguilar-Rodriguez_2016 | +0.022 |
| Bustamente_2010 | +0.038 |
| Payne_2016 | +0.059 |
| Arizaga_2000 | +0.075 |
| Molina-Freaner_2003 | +0.132 |
| McMullen_2009 | +0.164 |
| Molina-Freaner_2004 | +0.183 |
| Lemaitre_2014 | +0.332 |
| Morse_1983 | +0.526 |

Interpretations: the global study-equal median +0.01103 differs from its study-equal **mean -0.14677** because several negative deltas are large; both require sensitivity to endpoint denominator, source layout and paper selection. The seed/formed-fruit numerator is a conditional post-fruiting trait, not total seeds per initially assigned flower, and raw plant-level pairing and source variances are absent. No plant fitness mediation fractions, fixed-effect meta-analysis, or causal pollination guild mechanism follows.

## Consequences for the ALAN manuscript

1. **Withdraw universal seed-stage amplification** as a main claim. The strongest new result is *heterogeneous stage allocation*; exact stage matters in syntheses and within-study (not only between-study) controls are critical.
2. Preserve the distinct ALAN reproduction signs: Knop negative fruit; Lu night-only conditional seeds negative but open seeds unknown; Macgregor and Li positive configurations; Grognuz daylight carryover; Boom conditional seed predation. They still do not form a variance-harmonized pooled light-treatment estimate.
3. The original **9 directly manipulated reproductive ALAN papers / 3 skyglow gradient papers** are unchanged. Source B 113 diel-access studies are NOT A-light studies. **ALAN reproductive lnRR × independently measured same-population dark V** remains **k=0 confirmation-ready**, beta not estimated.
4. New natural next test: recover original *initial flower → formed fruit → viable seeds* data within multiple independent lighting experiments before claiming mediation; if no raw denominators, maintain a systematic evidence map and source-request programme.

### Method status

Post hoc original study reanalysis, not a new PRISMA-complete search. Paper bootstrap is only a descriptive resampling interval: original arm-level SE, within-study covariance, phylogenetic nonindependence, shared research-group selection and bias are not modeled. Same study, site, year, period and other treatment **does not guarantee** the two endpoints were measured on the identical plants or initially tagged flowers. The apparent stage difference is NOT a causal effect of the night pollinator guild and absolutely not a confirmed ALAN reproductive response.
