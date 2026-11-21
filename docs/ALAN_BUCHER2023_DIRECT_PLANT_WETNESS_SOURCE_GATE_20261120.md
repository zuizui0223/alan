# ALAN mechanism bridge: original plant biomass × moisture source support

2026-11-20: **SOURCE ACCESS STOP**, no plant×MOIST coefficient fitted, \`main\` unchanged.

## Research goal

Directly connect the macro finding that long-term ALAN-associated recording-decline candidates have higher wet-habitat affinity to an *independent* manipulated-light experiment with **multiple plant species**.

[Bucher et al. 2023, Philosophical Transactions B](https://doi.org/10.1098/rstb.2022.0358), [official Dryad dataset DOI 10.5061/dryad.hhmgqnknt](https://doi.org/10.5061/dryad.hhmgqnknt). The 2023 Ecotron model communities contained 16 sown species exposed to a gradient of nocturnal illuminance. The public repository lists \`ALAN_specis_biomass.csv\`, \`ALAN_biomass.csv\`, \`ALAN_traits.csv\`, and \`README.md\`; source descriptions say \`unit\` = Ecotron, \`light\` = lux, \`mass\` = species-specific dry weight, \`harvest\` = 1 or 2. Multiple harvests/species within an Ecotron are **not independent lighting replicates**.

This design could test whether independent European plant MOIST affinity predicts species-specific *experimental* biomass response to illumination, using Ecotron as treatment allocation unit and source-verified species name mapping. It cannot test true wetland drainage because all plants shared experimental soil.

## Source retrieval completed with a valid STOP

[Successful SHA/schema/access-audit CI #44657052979](https://github.com/zuizui0223/alan/actions/runs/44657052979) ran the [source-gated original data fetcher](../scripts/acquisition/alan_bucher2023_dryad_original_species_gate_20261120.py). All four official file IDs were checked at three source endpoints each. On 2026-11-20, official Dryad \`downloads/file_stream\` returned HTTP **403**, and the published \`api/v2/files/{id}/download\` endpoints returned HTTP **401**. No authentic CSV bytes could be retrieved or validated in this execution. Exact [access receipt](../results/all_species_forest/biotic_interactions/bucher2023_dryad_species_source_access_20261120.json).

Accordingly:
- No EIV exact species match has been made to original species biomass row labels.
- No source 12 EcoUnit (or alternate) sampling frame, actual intensity levels or per-species response was independently verified. The original paper's study design is not a substitute for source row QA.
- No nightlight × MOIST, direct plant stress, or independent wet-habitat mechanism coefficient has been estimated.
- No unverified source numbers or hand-transcribed published summaries are presented as new analysis.

## Once genuine source files are independently accessible

The frozen next-step specification is to verify original archive SHA, CSV columns, treatment assignment, Ecotron unit IDs, two harvests and species abbreviation-to-full-name mapping **from authors' supplied tables**. Match raw EIV exact species names, keep unknown traits unknown, and require wet and drier species coobserved in enough EcoUnits. Fit a species×light association on a plot-level experimental unit while retaining species × year/harvest and repeated unit structure; use *unit-level* randomization/bootstrap rather than species rows as independent units. Pre-specify light dose, outcome transform, null/no-response, diversity and individual plants as separate endpoints. Source missing means **STOP**.

This would be a candidate direct-plant physiological/competition mechanism, whereas Lu 2026 field experiments and Page 2021 network q×V explore biotic compensation. Neither experimental strand by itself proves the original wetness-biased long-term record changes.
