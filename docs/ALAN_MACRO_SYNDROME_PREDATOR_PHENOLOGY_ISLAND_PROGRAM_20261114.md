# ALAN — A macroecology of nocturnal pollination syndrome vulnerability, predator dual roles and phenological mismatch

**Source-backed status: 2026-11-14 JST.** Real archive-aggregate QA and predeclared future estimand gates passed: [GitHub Actions #43564470445](https://github.com/zuizui0223/alan/actions/runs/43564470445), [aggregate artifact #14101791804](https://github.com/zuizui0223/alan/actions/runs/43564470445/artifacts/14101791804). This is a macro-programme **with an executed source-coverage validation**, NOT a newly estimated causal global interaction coefficient. Prior ALAN positives and nulls are preserved unmodified; parent PR #63 and main untouched.

## Central ecological question

**Which pollination syndromes and combinations of plant reproductive dependence, animal community substitution, predator pressure, inherited island pollination channels, and flowering–pollinator synchrony make plant reproductive functions sensitive to artificial light at night?**

The appropriate unit is not "individual moth-visited plant species" but a **plant-functional guild or verified pollination syndrome × plant population × physical island/mainland locality × year × season × daily time**. Study/paper and plant phylogeny are independent statistical grouping layers, not ecological treatment units. Associations of a floral shape/color with light are useful preliminary observations but cannot prove true moth reliance. Records missing from GloBI/GBIF are UNKNOWN and not absence.

The core proposition is **channel-gated vulnerability**: night light will create reproductive losses when light-responsive effective nocturnal pollination becomes lost or unsynchronized AND the plant has insufficient alternative daytime/vertebrate/lizard/insect pollen delivery and autonomous seed assurance. Negative effects may be offset or reversed if predators reduce florivores, if animals exploiting moths also directly deliver adequate pollen, or if alternative pollinator guilds increase **effective conspecific pollen deposition**, not simply visit frequency.

## Four explicitly distinct mechanistic hypotheses, not automatically a four-way fitted regression

### H1 — ALAN × realized night dependence × reproductive assurance
- Measure same population/year seed or fruit response from night-only, day-only, open and valid closed treatments (plus hand pollen addition where its methods are actually comparable).
- Confirm **effective** nocturnal pollen transfer, and distinguish physiological self-compatibility from measured autonomous selfing and pollen limitation.
- Predict sensitivity to light may be strongest in self-incompatible / demonstrably nocturnal-dependent plants with weak alternate pollinator pollen delivery.
- A 2026 original field experiment [Lu et al., J Plant Ecol](https://doi.org/10.1093/jpe/rtag160) measured a 56% seed-set decline in nocturnal-dependent self-incompatible Pedicularis, but no detected response in a self-compatible diurnal/selfing Euphrasia. Only two species and one field experiment: NOT a global interaction.

### H2 — ALAN × predator pressure × sign of service
Two competing pathways require measurement at the same place and time:
1. **Negative:** adult moth interception, predator-induced foraging fear, reduced actual pollen transfer; predators may also intercept butterflies or beetles.
2. **Positive:** predators remove florivores and/or a nectar-feeding gecko, bird, bat or other animal actually transfers conspecific pollen and thereby **compensates** for pollinator losses. A flower visitor is not automatically an effective pollinator.
The same species may have both positive and negative potential roles with different prey stages and flower resource conditions; its **net** contribution must be estimated, not assigned by taxon label.

The full 24,453,954 GloBI interaction source records six bat species with source-restricted moth-family predation + database pollinates relations; ZERO gecko and centipede species meet that strict independent dual-role criterion. These same-species names **do not** establish a local simultaneous role, larval/adult prey discrimination or effective pollen transfer. The 2013–2025 GBIF same-island/animal-photo observations support 15,255 common cells, but neither capture nor pollen delivery.

Source-authenticated original field bridges now include:
- Wenzell 2024, **one publication**, six linked real 2019 Castilleja population-year guild/fruit comparisons (5 within 1 km, 1 provisional ~1.65 km), published guild fruit function already known to authors.
- Kays 2012, **one publication**, original 2010 Gamboa *Ochroma* flowers with fruit treatments and kinkajou feeding on the actual night-only fruiting flowers; 'complete exclusion' original meta source in fact admits bees through wire mesh (no autogamy-only inference). Not a night-light manipulation.
- Marques 2015, additional independent paper with bat/hummingbird source visitation and original day/night fruit/seed outcomes; independently matched visitor year remains unknown.

### H3 — ALAN × daily floral availability × seasonal phenological mismatch
**Two non-equivalent temporal scales**:

Daily overlap: proportion of reproductively available flowers at a local hour h (including opening/closing, nectar/scent activity) versus **effort-standardized effective pollen delivery** by actual local nocturnal/diurnal pollinator guild h. Use cyclic day/night observations, not just "night blooming" source labels.

Seasonal overlap: flowering onset/peak/end by day-of-year versus pollinator presence/activity distributions during the same locally observed breeding season, controlling temperature, daylength, climate and effort. A shift in photo upload dates or iNaturalist observation dates is **observer recording phenology**, not automatically biological first flowering or pollinator mismatch.

Example common-support overlap diagnostic, for normalized nonnegative profiles f and p on EXACT same time bins:
\[
O(f,p)=\sum_{t\in T}\min(f_t,p_t).
\]
Compute separately for lit/dark sites with identical time support; compare **within-year** seasonal/diel differences and actual fruit/seed outcomes. Effective pollen deposition weights can replace visitor count where measured. Since this index depends on sampling/effort and profile normalization, additionally retain absolute flower and pollen resource supply. Do not claim a mismatch solely from flowering shifts with no pollen vector records.

Original Swiss LED-field experiment [Grognuz et al. 2026, Biol Conserv 314, 111650](https://doi.org/10.1016/j.biocon.2025.111650) compared 14 lit and 14 control sites; 75% of 16 species showed seasonal schedule shifts and some plant opening/closing times changed, with cases >10 days. It does not independently measure same site/year pollen-vector mismatch or seed consequences. The global ALAN 2012–2025 photographic panel (332,576 source photo site-years) **is not** a first-flowering or effective pollinator phenology dataset and has no original official annual VIIRS for early-photo years.

### H4 — island × mainland × ALAN × *previously realized* channel shift
Island isolation does not force nocturnal moth shift in every community. Existing systems can differ among moth, bird, bat, lizard, introduced honeybee, fly or mixed generalized pollination due to actual pollinator supply, predator history, plant lineages, introduced species and island size/isolation.

Pre-existing pollinator channel shift is a **historical population-level state** confirmed with direct visitation/pollen or ancestral mainland/related phylogeny, NOT a trait assumption based on islands alone. Distinguish newly lit populations from those already reorganized before lighting; include ecological baseline selection. Keep island naturalness, size, latitude, mainland distance, source floristic quality, human/urban light, invasive vs native geckos, actual flower/bird/bat service separate.

Original photo-plant × same-year VIIRS island/mainland WLS: **10,698 site-years, 933 island/9,765 mainland, 78 geoblocks**, GloBI v0.11 two positive documented guilds vs moth-positive only **β−.1376 pp/loglight, CI [−2.048,+1.773], p=.886**. This is *source-documentation composition*, not actual shift. An earlier five-year legacy GloBI −1.9017 slope reverses to +.750 after v0.11 record update; not a biological nocturnal guild loss. 562-island historical PCNL only includes **52** photo plant species with source night/day reproductive evidence, so true shift×ALAN at species-site-year is not currently identified.

## Source evidence QC executed in this branch

Original signed aggregate evidence (not a cross-source synthetic 4-way regression):
- **562 islands**; 77,647 photo site-years; 2,012,925 accepted plant×site×year events; 14,330 recorded accepted plants.
- **52/14,330 = 0.363%** photo species with external real night/day reproductive evidence. Only 1,691/14,330 (**11.80%**) have measured source color and form; 257 white/tube-like, 13 also documented deep tube. The 52 and 257 strata are NOT known to share the same members; their intersection cannot be invented.
- 301 independently reproduced actual day/night experiments /118 papers /123 plants, observational earlier historical light association β+.003068 (95% CI −.02239 to +.02852, p=.813).
- 15,255 photo-only gecko/moth same-cell years; gecko share × historic light vs moth share β−1.386 (CI −5.457 to +2.685), unresolved.
- Full GloBI predator/pollinator same-name recorded relationships across guilds, but zero verified site/year predator-attack **and** pollen-delivery dual-role events under manipulated light from these source archives.
- Source compilation **p-value transcription bug**: results/alan_v011_moth_butterfly_flower_documentation_multi_channel_ALAN_20261112.json has \`frozen_positive_result_unmodified.p_nominal=0.4799746115\`, copied from updated **PL30**. The verified source original **PL64 nominal p=0.0439934091**, Holm p=.08799. Both original records and the historical source manuscript remain unchanged; this report identifies the provenance discrepancy.

**Independent new paper experiments supporting *parts*, not the full H1–H4 chain:**
- [2026 Tibetan experimental pollinator network](https://doi.org/10.1002/ece3.74368): 5 manipulated field blocks, 29 night-blooming plant taxa / 2,160 night visits, 40 day taxa / 3,381 daytime visits. Night Lepidoptera **−57.5%**, nocturnal Diptera **+179.6%**; day taxon composition changes even when total daytime visits do not. No plant-specific mature seed set in this network study.
- [2026 reproductive ALAN experiment](https://doi.org/10.1093/jpe/rtag160): 2 plants with contrasting reproductive assurance, same study group/area as network study; should not be counted as an independent global replication of regional ALAN effects.
- [2026 Swiss plant phenology experimental study](https://doi.org/10.1016/j.biocon.2025.111650): 14 lit /14 dark site pairs, measured daily and seasonal flower schedules; not paired pollinator synchrony.

## Interpretation and inferential gates

The executed real-source coverage QA status is:
- **H1: PARTIAL_EXPERIMENTAL_EVIDENCE_NO_GLOBAL_COMMON_SUPPORT**
- **H2: NOT_ESTIMABLE_DIRECT_PREDATION_AND_COMPENSATORY_POLLEN**
- **H3: NOT_ESTIMABLE_DIRECT_SEASONAL_OR_DIEL_MISMATCH**
- **H4: NOT_ESTIMABLE_EFFECT_OF_PRIOR_ISLAND_POLLINATOR_SHIFT**

Predeclared gate: at least **25 genuinely independent source papers and 20 independently documented plant taxa** in each syndrome group for inferential global cross-study comparisons, plus site/year exposure and relevant function measurements. This is a *minimum evidence gate*, not a power guarantee. A 562-island/photo opportunity frame and 24.45m source interaction claims cannot satisfy missing functional and temporal links by sheer row count.

## Executable next steps without changing ecological endpoints

1. **Experimental evidence atlas:** Extract the published plant×treatment×site×year×day/night outcome and *actual lux and spectrum* (not remote PCNL), with full methods for closed and hand access and inferred confidence limits. Expand field reports beyond Kendal/Nicholson and beyond just 2 plants.
2. **Functional channel map:** At same locality, separate plant nocturnal dependency, guild-specific effective pollen, confirmed autonomous selfing, actual predator attack and florivore removal. Keep predator's direct floral services independent of prey abundance until both observed.
3. **Seasonal/diurnal phenology matched time bins:** Obtain hourly flower opening times or daily flowering stage, pollinator captures/visit exposure by hour and calendar day, and fruit set. Infer overlap and stage shift under matched lit/dark, not from observation month of opportunistic photographs.
4. **Island–mainland matching:** Build island pollinator-shift baseline from source-confirmed history of effective pollinator guilds, restrict to climate/region and plant lineage common support, then test light interaction. No false control by calling undocumented mainland nocturnal alternatives biologically absent.
5. **Hierarchical analysis:** First fit ALAN × true-night-dependence with reproducible study-block effects, then ALAN×actual functional-redundancy and ALAN×phenology overlap where separately estimable; predator-negative and positive service paths only where observed. Attempt island×prior-shift×ALAN interaction ONLY when independent source support permits. Hold out whole source studies/island groups and control multiplicity; document all nulls.

The novelty is not re-reporting existing arthropod photographs or source categorical moth claims as proof. It is a **causally discriminating macroecological synthesis of nocturnal pollination vulnerability conditioned on channels, predator services and temporal mismatch** with transparent support gates and an executable original-source receipt.

## Reproducibility

Original source QA runner: \`scripts/analysis/alan_macro_syndrome_predator_phenology_source_gates_20261114.py\`, tests: \`tests/test_alan_macro_syndrome_predator_phenology_source_gates_20261114.py\`, predeclared source ledger: \`config/alan_night_syndrome_multichannel_macro_sources_20261114.json\`, CI: \`.github/workflows/alan-macro-syndrome-predator-phenology-source-gates-20261114.yml\`, SUCCESS #43564470445, aggregate artifact #14101791804. Source-only aggregate; no private plant/site/year coordinates or unverifiable ecological absence invented.
