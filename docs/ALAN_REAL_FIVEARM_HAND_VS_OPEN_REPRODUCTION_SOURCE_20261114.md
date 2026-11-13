# ALAN — true original five-arm HAND reproductive headroom audit

**Original-source GitHub Actions SUCCESS** on 2026-11-14 JST: [run #43529862353](https://github.com/zuizui0223/alan/actions/runs/43529862353), [aggregate-only artifact #14086129825](https://github.com/zuizui0223/alan/actions/runs/43529862353/artifacts/14087285335), [Draft PR #59](https://github.com/zuizui0223/alan/pull/59).

## Question and source

Does previously observed weak original four-arm DAY/NIGHT/OPEN/CLOSED reproductive nonadditivity survive comparison with a fifth original HAND pollination arm, an exploratory proxy for potential reproductive headroom?

Original published source: Kendall & Nicholson (2025), *Ecology Letters*, DOI [10.1111/ele.70036](https://doi.org/10.1111/ele.70036); 1,652 treatment rows, Git blob SHA 39a82d966afcdde670603add1c881ace8791ffe9. The source has direct day/night reproductive treatment outcomes, **not manipulative lighting**. Original study×plant×site×experiment-year×season×fruit-or-seed outcome units and source treatment means are exactly matched; missing HAND is never recoded as zero. True prior 1992–2021 historical PCNL satellite shards remain SHA/MD5 verified but are not reinterpreted as a light intervention.

The original author paper **already reports** flowering-time associations, day-night flexibility and ordinary pollinator dependence. This is a limited original-data reanalysis, not discovery of those published basic phenomena.

## Source authenticity and results

- Original field source: **301** matching DAY+NIGHT pairs from **118** papers / **123** plant species.
- Original full four arms: **153** same-unit pairs / **79** independent papers / **85** plant species; original four-arm nonadditivity study-median **I=−.05747065**, 49/79 negative, nominal p=.042165; historical parent two-endpoint Holm-adjusted p=.08433.
- Exactly **63** independent-metric fifth-arm pairs have a valid additional original HAND outcome, representing **32** independent source papers and **39** species, no zero or missing substitutes.
- Let H=original hand-treated reproduction and O=original unrestricted open reproduction; exploratory **headroom P=(H−O)/(H+O)**. HAND exceeds OPEN in **50/63** genuine matched outcomes. Across original paper-level medians **26/32** are positive, exact sign **p=.00053505** (nominal, posthoc). Source-paper P median **+.13377261**; original outcome P median **+.11054638**.
- In the same fifth-arm subset the actual four-arm nonadditivity I paper median **−.04780997** with **20/32** negative papers; exact sign **p=.215327**, not a supported independent pattern.
- Independently aggregated original plant taxa: **31/39** plant medians P>0, exact nominal exploratory sign **p=.00029408**. These 39 taxa share evolutionary history, and the data were not randomly selected from global plants.
- Study-median Spearman **rho(P,I)=+.1085044**, p=.554438; 3,000-study bootstrap descriptive 95% CI **[−.27972,+.46577]**; no resolved relationship.

The original source independent desk check spans 14 reported countries, 1993–2021 experimental start years, and 62/63 source pair records reported native plant status. This does not make the sample globally random or taxonomically representative.

## Key limits

Original HAND protocols are not uniformly standardized to a measured pollen-saturated maximum. Greater HAND than OPEN is a real published-treatment contrast, but cannot automatically be called a universal pollen-limitation estimate. The P and I ratios both share the same OPEN arm, creating algebraic correlation even without biology.

Outcome-selected P>0 (50 comparisons, 28 papers) and P≤0 (13 comparisons, 10 papers) descriptive subgroup summaries are **not independent confirmatory tests**: some papers contribute both; the latest script suppresses subgroup p-values. Fully measured five-arm source remains relatively small.

No causal source-verified link to artificial nightlight, local moth losses, reproductive insurance, autogamy or island floral change has been demonstrated. PR #51 authentic night reproductive share vs prior PCNL remains unresolved (p=.813); PR #58 prior PCNL vs four-arm I unresolved (p=.906). GloBI documented moth/butterfly positives do not identify pollinators at original experiment site/year (PR #57).

## Reproduce

Source-locked GitHub Actions #43529862353 SUCCESS. Code path: scripts/analysis/alan_original_diel_fivearm_hand_headroom_20261114.py. Tests: tests/test_alan_original_diel_fivearm_hand_headroom_20261114.py. Workflow: .github/workflows/alan-original301-diel-fivearm-hand-headroom-20261114.yml. Genuine aggregate artifact #14086129825. Original study×species×site/year rows and plant traits are not published in this fork. Keep this stacked PR Draft; all earlier PRs and main unchanged.
