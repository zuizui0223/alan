# ALAN — Gecko × moth × artificial nightlight: genuine global source and empirical evaluation

2026-11-10, stacked draft. Source class: GBIF **documented photographic observation composition** and source-corrected satellite nighttime light, not direct night predation or floral pollination.

## Ecological hypothesis

Artificial light may concentrate moths around luminaires, increasing encounter opportunity for geckos. If geckos exploit that subsidy, lighting may produce (a) higher gecko foraging activity at lamps, (b) apparent local prey concentration in some observation methods, (c) depletion or redistribution of moths away from flowers or dark refugia, and possibly (d) modified nocturnal pollen delivery. These steps are separate estimands; no source here directly measures gecko predation on moths, flower visitors, seed set or true density.

Existing field evidence: Japanese house geckos repeatedly return to lit microhabitats (Kobayashi & Mori 2020, https://doi.org/10.5358/hsj.39.184), introduced Hemidactylus frenatus in Australia uses brighter sites than native Gehyra dubia (Zozaya et al. 2015, https://doi.org/10.1111/aec.12287), and 2026 Cambodia surveys report 28% occupancy of illuminated advertisements by H. frenatus and higher presence at stronger illumination (https://doi.org/10.1111/jzo.70102). Some predominantly diurnal Phelsuma use lights for nocturnal activity (2021, https://doi.org/10.1111/aec.13012). Geckos are not universally light-attracted or exclusively insectivorous. Individual Gekko japonicus may exhibit weak light-site fidelity and light-preference learning results do not justify simple prey-driven causal claims (https://doi.org/10.1016/j.beproc.2022.104806).

## Source-frozen global data actually analyzed

- GBIF frozen original 2026-08 snapshot DOI https://doi.org/10.15468/dl.8yrbe7, genuine 11,634,053 Squamata original records (3.1014 GiB). Real source family taxonomy, year, coordinate, basis of record and precise GSHHG 2.3.7 polygon membership re-read from original Parquet. Original row-quality audit: https://github.com/zuizui0223/alan/actions/runs/42791325209. 2013–2025 human photographs on 20-km2+ GSHHG islands: **66,038 source photo observations of six named gecko families across 581 actual islands**, 23,024 gecko-bearing island×0.05deg grid×year cells; **35,361 Scincidae control lizard records across 448 islands**.
- Original genus-subset recording totals: Hemidactylus 19,269, Gehyra 1,749, Gekko 10,116, Tarentola 7,401, Phelsuma 13,238, Lepidodactylus 2,896, Cyrtodactylus 1,633. Names are field-source photo records, not standardized abundances or native/invasive assignments on each island. Pygopodidae excluded from main gecko family grouping; day geckos such as Phelsuma not automatically nocturnal.
- Independently verified source moth and butterfly families (Noctuidae/Geometridae vs Nymphalidae/Pieridae), historical real PCNL year-minus-one 2012–24: **62,688 source insect photo grid×year observations across 196 source islands**, from 22,713,063 real indexed 4-family records on 1,217 islands overall. Main source 2013–25 method defined in ALAN PR #37.
- Source-lagged co-support: 2015–2025 insect cohort **194 independent islands and 52,995 original photo site×years**; **97 islands with gecko-family photographs in some earlier source year**, 460 original within-island-year pregecko-photo-documented light-contrast groups, 533 groups without an earlier gecko photo. Scincidae-prior observation positive in 354 groups. Source year strictly earlier than insect record year, but because source starts in 2013, pre-2013 gecko presence cannot be evaluated.

## The actual first interaction result (not mechanistic proof)

Response Y is the **source GBIF photographic record fraction** of the two named moth families among the four original moth+butterfly families, NOT actual counts of moth individuals, feeding, effective pollen transfer or persistence.

Source-tested within-original-island × insect-photo-observation-year fixed effects, nominal true-island-cluster uncertainty, independent island-equal photo-site weighting, adjusting original 4-family source documentation volume and latitude. Main gecko-photo documentation × year-before-original-PCNL interaction: **β = +1.6922 percentage points**, nominal original island-cluster 95% CI [−0.4538,+3.8382], p=.1215. Thus NO clear differential light slope in gecko-photographed islands. The earlier-photo undocumented reference is NOT a biologically gecko-free control.

Scincidae prior-photo documentation × historical light interaction β +0.5215 [−1.7428,+2.7859], p=.650. Hemidactylus historical-source interaction β +1.6922 [−0.4622,+3.8466], p=.123. Phelsuma supports only nine source islands in the eligible source group, so no credible formal genus-specific slope test.

Latitude-adjusted and latitude+area-adjusted slope interactions for gecko photo documentation became nominally positive (β +3.213 p=.0123 and β +3.233 p=.0119, respectively), but the prespecified simpler model is null, country/latitude/observer effort confounding is substantial, and other weighting/era/source-exposure choices are null. These exploratory nominal alternatives MUST NOT be selected post hoc as proof of predation or beneficial pollination.

Most importantly, conditioning on ISLANDS WITH ANY PRIOR reptile photographs leaves **97 gecko-documented islands and only 7 gecko-undocumented islands**; this comparison FAILS the independent-source positivity gate. The usual island population in which there is no gecko GBIF photo includes many unsurveyed lizard sites. Lack of a photo is NOT evidence of gecko absence.

[Original observation-source and historical-light model SUCCESS](https://github.com/zuizui0223/alan/actions/runs/42792547085), signed numerical details ../results/alan_global_gecko_x_priorALAN_moth_photo_observation_20261110_verified.json.

## Actual original 13-year 5-km dark-land geometry interaction

To test whether gecko prior source documentation modifies an association with available dark terrestrial habitat separately from brightness, authentic 2012–24 PCNL 5 km island-land-minus-lake dark percentage and four-connected patch index were joined to exactly the same original insect photo site/year source. The original 62,688 source photo sites had 62,685 valid historical land pixel neighborhoods; two source thresholds 0 and 3, not independent replicated studies.

At threshold PCNL=0, adding dark-land/connectedness support and 2015+ source records retained **161 actual independent islands / 45,938 observed photo source cell×years**. Light×historical gecko documentation β +2.114 [−0.774,+5.003], p=.150; **dark-land fraction×gecko documentation β +10.054 [−5.772,+25.880], p=.211**. Neither resolves a modifier. Scincidae light/dark interactions also have CIs spanning zero. At PCNL=3 the same real source gives 166 islands / 47,781 cells, light×gecko β +2.272 [−0.628,+5.172], p=.124; dark-land×gecko β +7.805 [−7.782,+23.391], p=.324. Dark patch continuity independent main slope is also unresolved.

These are PHOTO-record associations only; adding gecko/land/ALAN terms and changing eligible source support changes the definition of the main dark fraction slope. The earlier source-only beta dark land +12.10 from another cohort must not be read as proof of moth benefits or changed by gecko predation. Within-island-year cluster and 10km non-overlap source selection code is included for geographic robustness, not independent population replicates.

[Actual 13-year gecko×darkland source model SUCCESS](https://github.com/zuizui0223/alan/actions/runs/42793170640), signed numeric summary ../results/alan_global_gecko_x_13yrPCNL_darkarea_moth_photorecord_20261110_verified.json.

## Why this result cannot establish the proposed mechanism

1. ORIGINAL GECKO SOURCE: GBIF human gecko photos are neither observed local live-gecko abundance, moth consumption, nor verified predator occupancy. A zero gecko record is unknown detection, not absence.
2. ORIGINAL MOTH SOURCE: moth versus butterfly photographed species/family is a measurement/observation composition that varies with observers, documentation, plantation/horticulture, urban lighting, camera use, and day/night photography availability. It does not count true moth abundance.
3. TEMPORAL ORDER: Source documentation gecko photographed at an earlier time is not proof the gecko was present when a specific moth specimen was observed, nor does it establish a causal light-mediated predator response.
4. LANDSCAPE SPATIAL RESOLUTION: original 0.05-degree photo cells and ~1km historical PCNL cannot identify individual 0–100 m streetlight attraction, direct gecko attack near lights or dark microhabitat at flower height.
5. EVOLUTION/INVASION: Hemidactylus frenatus is invasive on some islands, native in others, but original source does not contain independently validated island-specific status. Phelsuma day geckos can shift into night activity; do not pool all gecko taxa as obligatorily nocturnal.
6. UNMEASURED OUTCOMES: no species-by-prey-stomach contents, camera predation counts, night-flower visitation, pollen deposition, seed set or pollination limitation. Both positive and negative source associations are therefore noncausal.

## Scientific judgment

The source and source-observation cross-taxon geography are globally tractable (**581 islands with gecko photographic records; 97 documented gecko islands in the current 196-island moth ALAN contrast**), and nightlight-dependent predation is a plausible mechanistic question supported by independent local field literature. But in the new macro GBIF-photo source analysis, the gecko×light and gecko×dark-land interaction coefficients are not stable or resolved under the source main/observer-coverage gates. A claim of gecko-mediated suppression of moth pollination is NOT currently supported.

Next evidence gate: validate island-specific native/introduced gecko identity from an independent authoritative distribution source, and independently measured *lamp-proximate* time-stamped gecko attack/foraging and moth flower visits or light-independent insect monitoring, on islands with and without earlier documented geckos. Treat any subsequent source-only interaction as observational until this is done.