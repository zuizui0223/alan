# ALAN | Moth aggregation, light-trap bias and effective pollination are different endpoints

2026-11-10 — true data-source status and testable causal alternatives.

## Hypothesis

Artificial light can **concentrate** moths near lamps even when **ground-level moth activity and pollen transport decline**. Conversely, special plant–moth mutualisms can respond positively to nightlight. Neither implies that flower-photo assemblages with more plant taxa documented in GloBI as moth-visited gained effective pollination.

Published matched light/dark results (Macgregor et al., Global Change Biology 2016, DOI 10.1111/gcb.13371): flight activity at the height of lights **+70%**, ground-level moth abundance **approximately half**, richness >25% lower, and overall pollen transport lower at lit sites. Published yucca–moth field work (Wilson et al. 2021, DOI 10.1093/icb/icab010): local adult moth trapping and plant fruit/offspring response can be positive under some forms of artificial nightlight. Knop et al. 2017 (DOI 10.1038/nature23288) found **-62% nocturnal flower visits** and -13% fruit set at experimentally lit plots.

## Test four separate quantities

1. **Aggregation:** moth counts or flight altitude near a lamp; positive phototaxis may increase these.
2. **Functional foraging:** moth visits to flowers, pollen carriage/deposition, or nectar-feeding; these may decrease, increase or remain unchanged, depending on species, environment and light architecture.
3. **Demography:** egg laying, larval recruitment, adult longevity or survival; local trap numbers do not indicate these.
4. **Plant population response:** flower/seed/fruit number and later recruitment. The ALAN global +0.449 photographed GloBI-claim fraction is in a completely different source class.

Never treat the first quantity as proof of the other three; never label absence of a GloBI source claim a known non-moth plant.

## Independent source promising for original quantitative re-analysis

Van de Schoot, Wesselingh & Van Dyck 2025, Environmental Pollution 372:126068 (DOI 10.1016/j.envpol.2025.126068), Dryad original DOI 10.5061/dryad.6djh9w1d2. Public metadata describes true moth captures using **light and bait traps**, distance to nearest streetlight in metres (NearestSL), number of streetlights within 400 m (SL400), skyglow, vegetation, garden area and built-up land. It documents individual-sample and pooled-per-trap Excel tabs. Comparing standardized within-method abundances and a light-by-trap-method interaction can diagnose **sampling-method sensitivity** without equating bait traps with an infallible density census.

## Actual GitHub data acquisition result, not a biological estimate

[Dryad native original-data preflight Actions #42764795919](https://github.com/zuizui0223/alan/actions/runs/42764795919): **technical source audit succeeded**. The official Dryad dataset JSON endpoint served authentic metadata (published April 1, 2025; versionNumber 4). However, the official public **full data ZIP download returned HTTP 401 Unauthorized**. No individual row, original Excel table or new insect coefficient was obtained, and no output was fabricated. We do not bypass the access control.

The executed source-audit code is at scripts/acquisition/alan_dryad_moth_lightbait_preflight_20261110.py, output artifact 13739132165. The original dataset can be re-evaluated only when its actual public bytes are accessible through an authorized official endpoint or supplied directly by the original data owner.

## Conditional analysis plan, when source is genuinely available

Require the actual individual sampling events from the original Excel sheet, with original location, season/night and trap type. Check number of spatially distinct sites and contemporaneous observations in the same study. Report site × trap-type effort before any fit, then independently model observed counts for light traps and bait traps under distance to streetlight, ALAN skyglow, urbanization, vegetation, and date/method effort; compare cross-method light slopes as an **observational detectability test**. Two trap types differ in capture efficiency and species selectivity, so method interaction alone does not identify attraction or abundance.

The proposed 0-20-50-100m lamplight edge hypothesis is not tested by original 0.05-degree photographed plant cells or ~1-km PCNL satellite products. A measured fine-scale streetlight response is distinct from within-island landscape dark-refuge fragmentation, whose source-verified separate macro pilot is [PR #32](https://github.com/zuizui0223/alan/pull/32).

**Scientific status:** local positive moth attraction is plausible and observed in literature; neither positive effective nocturnal pollen delivery nor causal explanation of the previous positive photographed plant guild fraction has been established in ALAN.
