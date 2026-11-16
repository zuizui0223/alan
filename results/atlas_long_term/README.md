# Great Britain plant distribution changes and satellite night-light change

Exploratory independent regional extension of the frozen five-year analysis. Retain the original −1.902 and the 64-species pollen limitation result in their existing roles.

## Sources and estimands

Plant Atlas 2020 (Stroh, Walker, Humphrey, Pescott & Burkmar, 2023) provides [10-km period distributions](https://doi.org/10.5281/zenodo.10998003) and [effort-adjusted British national trend ensembles](https://doi.org/10.5281/zenodo.11108832). Annual harmonized DMSP-equivalent satellite night-light DN from [Li et al., version 10](https://doi.org/10.6084/m9.figshare.9828827.v10). Species-rank atlas entries without qualifiers or explicitly sensu stricto are retained. Other ranks, aggregates and ambiguous broad entities are excluded. Plant reproductive traits come from PRs 84 and 85, matched by scientific name, without assigning aggregate traits or filling missing information.

National decline is defined here as at least 95/100 publisher ensemble gradients in the moderate or strong decline bins (slope below −0.001). The 95/100 cutoff is an analyst choice. Apply the original minimum 16 British squares and remove records flagged for higher-aggregate substitution. Short estimates cover 1987–2019 and long estimates 1930–2019. These overlapping windows measure net relative distribution changes, not uninterrupted annual losses. They share records with the spatial light analysis.

For the spatial test, select OSGB hectads with at least 50 retained species recorded in each endpoint survey period. Each species starts from squares in which it was actually recorded during 1987–99. The response is whether it was recorded again during 2010–19. Colonization of new cells is outside this estimand. Average 25 regularly spaced light samples per square annually, then contrast 1992–99 and 2010–19 period means. A square can include sea. DN is a harmonized satellite index with saturation and a sensor transition; local lamp spectrum, lux and light before 1992 are unavailable.

Species linear probability models include initial light, log early/late richness excluding the focal species (recording proxies), and 100-km OS-square fixed effects; use finite-sample-corrected block-cluster standard errors and t reference. Require at least 50 baseline cells, 10 blocks, 10 redetections and 10 non-redetections. BH across all 1,492 estimable species. Thresholds are exploratory. Report observed-species interquartile differences alongside standardized DN coefficients because light changes span very different ranges. Richness proxies cannot separate observer turnover from genuine community changes; land-use, climate and exact observation effort are unresolved. Both light and persistence can covary with those factors.

## Outputs

- [Species-level joint evidence](species_evidence.csv) retains all 3,016 source species entries, including unknown values. 2,092 had usable short-period national estimates, 611 met the decline definition and 340 also met it for the long window. Of 553 declining species with a light coefficient, 259 point estimates were negative. None passed BH q<0.05 in that direction. The three FDR-supported negative associations in the entire light screen did not meet the national decline criterion.
- [Comparison including the reference species](trait_comparison.csv). The reference group comprises all other species with eligible short-period trends. Differences are descriptive; study coverage is highly selective. Counts of pollen transport and pollen limitation are restricted to previous shared-cohort evidence. Nine of the 611 have documented moth transport in that ledger, eight have a pollen-limitation estimate. Three have directly matched day/night experiments, all Caryophyllaceae. This reflects experimental coverage and cannot characterize the 608 unknowns.
- [Per-species coefficients and support](atlas_light.csv); [period light and recorded plant richness by public 10-km grid square](grid_support.csv). The 100-km block is not an island/mainland classification.

## Replay

Extract the two cited atlas archives under CACHE. Use the harmonized annual raster files for 1992–99 and 2010–19 as restored by alan_decadal_light.py. PR85 is the source checkout of 8c0c610. LEDGER is its checksum-verified direct_species_trait_ledger.csv.gz. OUT is results/atlas_long_term.

    python scripts/analysis/alan_atlas_rerecording.py CACHE LIGHT_DIR OUT
    python scripts/analysis/alan_atlas_declines.py CACHE OUT PR85 LEDGER
