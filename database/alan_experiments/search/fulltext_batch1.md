# Full-text eligibility batch 1 — 2026-11-20

Seventeen cached Europe PMC full-text XML records were individually screened using Methods/design, Results, and data availability or supplement statements. Decisions were not based on titles/abstracts alone. `fulltext_batch1.csv` records study identity, treatment unit, sample-size qualifications, source sections, and possible programme overlap. Empty replication cells mean unresolved or inapplicable, never zero.

| Decision | Records |
|---|---:|
| Eligible design, extraction pending | 4 |
| Hold: independence or scope | 8 |
| Separate indoor behavior | 2 |
| Separate observational | 1 |
| Separate spectrum | 1 |
| Reference/infrastructure | 1 |

**These are eligibility decisions, not 4 completed inclusions. No effect means, SDs or original-scale sampling variances were extracted in this batch.** Original datasets cited below were identified in primary full texts; availability statements do not by themselves prove usable raw independent replicates.

## Highest-priority interaction sources

1. **PMC4375371, doi:10.1098/rstb.2014.0131.** Outdoor grassland experiment: 54 mesocosms, white/amber/dark crossed with plant-only, herbivore and predator-containing communities, six replicates per cell. Supplementary datasets are explicitly provided. This is directly useful for distinguishing plant-resource changes from top-down predator hypotheses. The Results do not support a detectable added predator effect in the selected models; preserve that counterevidence. Later Penryn grassland publications require experiment/time overlap checks.
2. **PMC7368033, doi:10.1038/s41598-020-68667-y.** Four illuminated and four dark meadows with *Silene latifolia*, *Epilobium angustifolium* and the pollinating seed predator *Hadena bicruris*. Illuminated and adjacent-dark sectors share meadows. Seed production per capsule includes null treatment results. Original data are author-request, not a confirmed open raw-data repository. It is especially relevant to the dual beneficial/harmful function of an interaction partner. Do not mix its field phytometers with its separate controlled flight assays.
3. **PMC12205020, doi:10.1007/s00442-025-05750-z.** Outdoor earthworm experiment 2: six mesocosms per light × earthworm-presence combination; litter removal, respiration, survival and surface activity. Worm replenishment on day 16 and omitted nights matter. Final litter-removal and some respiration contrasts include nulls. Data are author-request; figures/supplement require assessment. The indoor repeated-exposure experiment 1 is not an independent field experiment.
4. **PMC7821510, doi:10.1186/s40462-020-00238-2.** Experimentally lit forest edges with red/white/dark conditions and bat flight trajectories. Original data/code: [Zenodo 4036279](https://doi.org/10.5281/zenodo.4036279). Four sites and sequential site-night sampling need row-ID reconciliation; the article's four sampled observations per condition must not automatically become four independent rows. Bat position or activity is not measured attack rate or prey removal.

## Mechanistically relevant studies that cannot enter the same lit/dark synthesis

Follow-up 2026-11-21: the [bat raw-data support audit](../bat_refuge2021_sources/README.md) recovered and reconciled 28,646 retained positions and 1,855 trajectories. The 12 recording nights comprise 10 site–treatment cells (control 4, red 3, white 3), with repeated red/white observations at one site. Reported coefficient inconsistencies are retained; no predation effect is inferred or added to the pooled catalogue.

- **Reflection/polarization, PMC6346644:** original [Dryad 7m5cv66](https://doi.org/10.5061/dryad.7m5cv66). Five river sites, three light spectra × two polarizing surfaces; all illuminated. Strong relevance to the water-reflection mechanism, but no dark reference. Keep in the separate spectrum/mechanism evidence layer.
- **Parasitoid retrofit, PMC11842737:** [Figshare 26206229](https://doi.org/10.6084/m9.figshare.26206229). Three road-section impact/control comparisons, all lit, combining shielding and spectrum changes. This is not strictly a wavelength-only study, so it is held as a mitigation-scope case rather than misclassified as a valid ALAN/dark experiment. Parasitoid capture is not measured parasitism.
- **Lakeside citizen-science attraction, PMC13594263:** [Figshare 29436527](https://doi.org/10.6084/m9.figshare.29436527). Three HPS and three shielded LED lamps in a nonrandom array, no dark arm. Repeated nights/time intervals do not increase treatment replication.
- **Mayfly beacons, PMC5717677:** two river-site crossover demonstrations with shore lamps remaining on. Adding overwater lights reduces attraction to existing shore lights, but neither arm is a naturally dark ecosystem. Repeated photographs cannot supply independent treatment n; mortality is inferred from attraction, not directly counted.
- **Oyster annual series, PMC12627416:** excellent high-frequency raw data on [PANGAEA](https://doi.org/10.1594/PANGAEA.983298), but one illuminated table and one control table. Oysters/time steps cannot replace independent illumination replication.

## Scope safeguards

The bird observational record (PMC11590875) was separated because ALAN was not manipulated, not because it studies birds. The other great-tit study (PMC4927185) uses experimental lighting but only one transect per color at a single site; logger-exposure data do not establish replicated ALAN population effects.

Indoor ecological aquarium growth/feeding (PMC10571867) is flagged for scope adjudication, not silently collapsed into the behavior-only exclusion. Indoor spider neural anatomy (PMC11371435) also needs endpoint/scope adjudication. Laboratory crayfish light preference (PMC12508872) and mosquito individual locomotor/metabolic assays (PMC9865375) are separated from the outdoor ecological synthesis. Field isopod manipulation (PMC13071524) exists, but lacks a replicated dark comparison; it is not mislabeled purely observational. LakeLab engineering (PMC8648721) establishes an illumination facility without biological effect estimates and is retained as infrastructure/reference context.

The broad question is therefore feasible as several explicitly different evidence types: replicated lit/dark ecological responses, illumination-mitigation experiments, mechanistic reflection experiments, and single-installation time series. Combining their effects as though they share a dark control or independent treatment unit would inflate evidence.
