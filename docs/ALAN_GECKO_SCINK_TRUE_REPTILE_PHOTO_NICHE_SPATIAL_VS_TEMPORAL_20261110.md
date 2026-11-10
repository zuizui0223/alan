# ALAN — Artificial light and gecko–skink photographic niche separation on islands

**Date:** 2026-11-10. Independent actual original GBIF/Prior-PCNL source research. Stack PR #39 on #38. SOURCE OBSERVATION COMPOSITION ≠ predator behavior or moth pollination.

## Scientific question and why the earlier model failed

Previous PR #38 compared gecko photo-documentation on 97 islands versus no earlier gecko photo on others, but when requiring actual earlier reptile photo opportunity, only **7 independent comparison islands** remained. A missing gecko photo cannot identify biological absence, and the main island-level gecko × light interaction was unresolved (β=+1.69, 95% CI [−0.45,+3.84], p=.122). Darkland × gecko-source interactions were also unresolved. This was a genuine source/detection-support failure, not a biological demonstration of no predation.

Here we instead compare **actual gecko-family photo records versus actually documented Scincidae (skink) photo records in the SAME island×survey year**. Both are genuine geolocated source categories, though both can have different seasonal, habitat and photographer visibility. The target is differential relative source documentation, not biological occupancy or abundance.

## Frozen original source and annual night light

- Frozen GBIF original 2026-08 DWCA preview DOI https://doi.org/10.15468/dl.8yrbe7, 11,634,053 original Squamata Parquet rows from 3 actual native files (approximately 3.10 GiB compressed), source-authenticated in PR #38. Filter original 2013–2025 georeferenced HUMAN_OBSERVATION PRESENT records using original GSHHG 2.3.7 20km²+ island polygons before aggregation to original 0.05-degree photo-cell × island × calendar year.
- **66,038 actual gecko-family photo source records, 35,361 actual Scincidae photos**, across **686 GSHHG source islands** in at least one of the two families, **34,056 source photo grid×year rows**. Original SHA256 source family grid CSV: c2e0c823ee9fcf61182b920724e488d9130b1d1ca9d53266511c2da8873983c9.
- Original corrected Li–Cao 1992–2024 v3 PCNL annual source DOI https://doi.org/10.5281/zenodo.17013414. Sampled each genuine reptile source grid *center* only when inside the SAME GSHHG island land, using 2012–2024 **observation-year-minus-one** genuine annually MD5-authenticated satellite raster. No 2024 night light assigned to earlier gecko photos, no ocean pixel made dark. All 13 years and four SHA-audited source workers succeeded. **28,814 true photo-cell years on 583 islands** have numeric preceding-year satellite source. Source run https://github.com/zuizui0223/alan/actions/runs/42796502910.
- True same-island×source-year contrast requires ≥2 separately located 0.05-degree photo cells, ≥2 genuine original satellite radiance values with log1p(light) span ≥0.02, and photo sources for BOTH gecko families and Scincidae in the island/year: **153 actual independent islands, 808 island×year blocks, 21,841 genuine source photo cell-years**, 44,227 gecko photo records and 24,768 skink photo records within these selected events. Stronger ≥5 original records for EACH family group: **80 islands/396 island×year blocks/19,308 photo-cell years**.

## Positive source spatial relationship, with cautions

Primary postexploratory fixed-effect model: original source gecko-photo percentage among the two true photo-guild classes ~ log1p(previous-year PCNL original pixel) + latitude of photo-grid center + log1p(original gecko+skink source event count), fixed intercept per **real island×survey year**, island-equal weighting, original island-cluster nominal 95% CI.

**β = +3.2669 photographic percentage points** per log1p(previous-year PCNL), CI [0.9115,5.6224], p=.00687 (153 independent islands/21,841 actual cells). In a stronger both-classes ≥5 original photo-events per island/year cohort, β +3.4674 [1.6253,5.3095], p=.000339 (80 independent islands). 10km-separated real source site centers (non-overlap of 5km geographic circles) yield β +3.5737 [0.1643,6.9831] (121 islands) and β +3.8481 [0.5388,7.1574] (122 islands). 2017+ β +3.6122; 2020+ β +3.0365; both nominally positive. Excluding the source most-photographed island, β +3.2596 [0.8798,5.6395].

However **photo sampling choices matter**: equal individual *source photo-cell* weights β +1.9101 [−0.6311,+4.4512] p=.140; at least 3 gecko/skink photo events per cell β +1.1853 [−1.3662,+3.7367] p=.358; adding area×light β +2.7183 [−0.3549,+5.7916] p=.083. The original source photo-event *volume* itself associates with historical PCNL: β +0.04670 in log(1+gecko+skink original photo count) [0.01935,0.07405], p=.00094. Thus an urban/photography/detection mechanism is plausible; the positive proportional relation is not automatically an ecological preference or trophic interaction.

All full 13 predeclared variants independently logged and archived in [signed source-result JSON](../results/alan_gecko_scink_spatial_positive_temporal_unresolved_20261110_verified.json), [original source-fit Actions](https://github.com/zuizui0223/alan/actions/runs/42797211231).

## Strongest counter-check: repeat exact original 0.05-degree grid across years

Within the stricter actual island/year source frame, retain only same GSHHG island×0.05deg site whose source has ≥2 distinct years, ≥0.02 within-site log1p preyear PCNL span, and at least one documented gecko **and** skink source across its recorded years. This leaves **91 original independent islands, 1,385 source photo-grid sites and 6,608 actual site-year events** (site repeated median 4 years). Regress gecko photo composition on historical PCNL and total source photo events with weighted **grid-cell fixed intercept + calendar year fixed intercept** and original island-cluster uncertainty, not merely island×year spatial regression.

**Island-equal β = −1.053** [−11.801,+9.695], p=.846. Source photo-cell/year-equal β +0.204 [−7.407,+7.815], p=.958. No resolved **within-same-site temporal gecko-share change** with source annual light. This is a much stronger spatial/temporal distinction than simply reusing one large source sample with repeated p values.

No strong inference is possible that turning on lights creates gecko predator attraction over time. The observation-time and land-use confounding remain; source annual satellite resolution is ~1km and photo cell label is ~0.05° (~5km), so individual 0–100m lamp zones cannot be resolved.

[Actual original repeated site/year Actions](https://github.com/zuizui0223/alan/actions/runs/42798154950), full signed numeric output linked above.

## What the result actually does and does NOT show

**Positive descriptive finding:** on 153 independently surveyed GSHHG islands with both reptile family photo classes represented in the SAME survey year, the spatial proportion of documented gecko photos tends to be higher at satellite-brighter source grid centres, with sensitivity to weighting and per-cell photo effort. This source association is compatible with artificially lit foraging opportunities, but also with housing habitat, built-up land, observers, photo camera/night capture, historical species turnover and gecko/skink differing diel activity.

**No temporal proof:** on 91 original sites repeatedly photographed across years, no resolved same-place annual association survives fixed site+year comparison. No live gecko density estimate, no documented gecko diet/moth capture, no light-mediated predator attack rate, no moth flower visits, pollen deposition or seed/fruit reproductive service.

Independent field literature:
- Kobayashi & Mori (2020), *Gekko japonicus* individual repeated visits to illuminated microhabitats, DOI https://doi.org/10.5358/hsj.39.184 (behavioral light site fidelity; not universal light preference).
- Kobayashi et al. (2023) laboratory *Gekko japonicus* work found no strong intrinsic light preference nor reinforcement by cricket-foraging near lamps, DOI https://doi.org/10.1016/j.beproc.2022.104806. This supports testing conditional rather than universally assumed attraction.
- 2026 Cambodian *Hemidactylus frenatus* lit-advertising-banner observational study reports local brightness-associated gecko occurrence and study-level site-surface-area effects, DOI https://doi.org/10.1111/jzo.70102.
- Predominantly diurnal *Phelsuma* spp. can engage in nocturnal activity under artificial lamps, DOI https://doi.org/10.1111/aec.13012, so Phelsuma is not a universal strict daytime negative control.
- Australia national wildlife light-pollution guidance notes invasive *H. frenatus* may exploit bright prey-rich spaces avoided by *Gehyra dubia*, but that specific competitive replacement requires native/introduced island-by-species validation and independent site-level behavior.

**Next measurement gates, not existing empirical results:**
- Verify actual eventTime and source timezone / diel observation opportunity; field presence of an eventTime column alone cannot establish reliable local solar time.
- Independent built-up cover, tree cover and observer/trap effort controls at matched source photo grid. Estimate condition on prior photos of both taxa, avoid gecko unobserved = absent.
- Original-night streetlamp-gecko-foraging videos, nearest light distance, moth prey capture and flower pollen deposition if the intended mechanism is gecko-mediated pollination disruption.
- Check if the gecko spatial association remains in previously moth-rich (independent observed insect) island-year source after source observer/detection controls. That would still be observational, not proof of consumption.

**Manuscript-grade phrasing:** "Across 153 GSHHG islands, records of geckos relative to skinks were spatially enriched under greater historical satellite nightlight in an observational source analysis, but the association was absent in repeated-site annual fixed-effect comparisons, leaving light-induced gecko foraging and moth-mediated pollination impacts unresolved."
