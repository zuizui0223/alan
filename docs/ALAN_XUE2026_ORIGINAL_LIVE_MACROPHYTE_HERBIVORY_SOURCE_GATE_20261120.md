# Xue et al. (2026): live-macrophyte herbivory as a candidate ALAN interaction mechanism

2026-11-20 • Source metadata authenticated, original biomass data **not** recovered. Draft source gate, \`main\` unchanged.

## Biological result independently documented in the original article

[Xue et al. (2026), Ecology 107:e70458](https://doi.org/10.1002/ecy.70458) manipulated artificial light at night (ALAN) and herbivore access around live macrophyte communities. In a 116-day factorial 2-light × 3-snail treatment mesocosm experiment, invasive apple snail *Pomacea canaliculata* reduced total native aquatic plant biomass by **48.0%** under no ALAN vs **87.2%** under ALAN (relative to same-light no-snail community baselines). The original paper's measured taxa are *Myriophyllum aquaticum* (nonnative in China) and three aquatic natives *Vallisneria natans*, *Hydrilla verticillata*, and *Myriophyllum spicatum*. Invasive *M. aquaticum* itself was largely avoided by the invasive snail; light enhanced invasion through differential growth and herbivory.

The implied fractions of native biomass remaining are 0.520 and 0.128, a ratio 0.128/0.520 = **0.24615** and a suppression contrast 87.2−48.0 = **+39.2 percentage points**. These values are **simple arithmetic of published author estimates**. They are NOT a newly fitted light×snail interaction coefficient, a standard error, or an independent source 57-replicate effect.

Experiment 1 original source design: **six physical outdoor lamp zones, three lit and three dark**, four source species/80 pots (77 after mortality). Source pots are **nested within assigned lamp zones**, not 80 independently assigned light installations. Experiment 2: source plant preference buckets, five replicates per listed snail×plant×light combination; physical light fixture hierarchy must be recovered before a field-independent-treatment claim. Experiment 3: 60 original mesocosms/57 final biomass, two lighting treatments × three snail treatments, ten nominal bucket replicates per cell. Original source article does not establish 60 independent lamp/illumination zones or actual allocation mapping.

## Original deposited ScienceDB data – API VERIFIED, FILE ACCESS STILL STOP

The author's Data Availability statement identifies **Science DataBank DOI [10.57760/sciencedb.35715](https://doi.org/10.57760/sciencedb.35715)**. DOI currently resolves to original dataset ID \`b9b02a9ec07c456b8584dd5fc732b7da\`.

We queried the official ScienceDB public API \`GET https://www.scidb.cn/api/sdb-openapi-service/json?doi=10.57760%2Fsciencedb.35715\`, whose route is documented by ScienceDB's OpenAPI description (\`https://www.scidb.cn/open-api/v2/api-docs\`). [Official GitHub Action #44661696531 succeeded](https://github.com/zuizui0223/alan/actions/runs/44661696531): HTTP 200, \`application/json\`, 3,870 bytes, SHA256 \`6e08a2cd77a5e663a98cdbbdce8aba91d73b273db0fe2a4d0925d4638686cbf3\`, with the exact published DOI present.

However, the DataCite/Schema.org-style JSON response does not list raw file assets or trusted original download URLs (no \`distribution\`/file manifests). **Zero genuine original 57-mesocosm rows, plant×treatment variance estimates, or light zone identifiers were acquired.** The original response was not mistaken for CSV data. [Source-original audit script](../scripts/acquisition/alan_xue2026_original_sciencedb_source_gate.py), [verified result](../results/all_species_forest/biotic_interactions/xue2026_sciencedb_liveplant_original_source_gate_20261120.json), [source tests](../tests/test_alan_xue2026_sciencedb_source_gate.py).

The original empirical **mechanism** result of light×selective grazing is persuasive within this focal aquatic community, but the manuscript's macro result that humid-habitat affinity is enriched among ALAN-associated recording decliners is a different study population. All four plants in Xue are aquatic; this trial cannot estimate whether the response gradient runs from dry to wet terrestrial plants, compare pollinator dependence or show ALAN-related reproductive failure.

## Quantitative next-step if authors' original raw file becomes available

1. Verify the DOI-matched raw dataset manifest/filenames/byte hashes; inspect original species, bucket, light-zone and snail-treatment keys before fitting anything.
2. For Experiment 3 primary response, estimate log biomass ratio of *native community in invasive-snail treatment* versus *same-light no-snail* separately under ALAN and no ALAN, then **difference in log ratios**, including covariance due to the physical light assignment. Treat the original **independent illumination allocation units** as the highest-level cluster, not 57 independent lights unless source layout confirms it.
3. Compare original plant-specific nonnegative biomass outcomes without filling unreported plant×bucket rows with artificial zeros. Distinguish native plant mortality and small biomass, herbivory intensity, M. aquaticum tolerance, and exposure-dependent invader growth.
4. Replicate in additional **independent** live aquatic and wet grassland light×herbivore studies, obtaining wetness trait variation in the **same experimental source**. This study on its own has no wet-vs-dry moderator or cross-study mean.

## Paper positioning

The original wet-species macro anomaly remains Results 1 (species composition, not ALAN causal decline). Xue's published selective herbivory is **candidate mechanism and external experimental evidence**. Without the raw data no independently reproduced interaction coefficient belongs in our manuscript Results. If raw source access becomes available, a new explicit effect size and correct uncertainty may be added; this cannot be described as the identified cause of the macro phenomenon without a wetness contrast.
