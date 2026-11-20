# Freshwater sediment communities: primary review, no independent light replication

[Primary Methods and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC4375370/?pdf=1), DOI [10.1098/rstb.2014.0130](https://doi.org/10.1098/rstb.2014.0130). Full HTML retrieved and cached externally as `microbes2015/article.html`; the bare web-tool endpoint failed but the alternate request returned the article.

The Westhavelland field comparison has one lit and one unlit ditch section, sampled before illumination and five months later. Three sediment samples per section are subsamples. The incubation places cores from both origins in one chamber: upper dark shelf and lower illuminated shelf. Temperature/photoperiod change sequentially from summer to winter. Field community composition showed a seasonal change without a significant overall illumination effect; laboratory oxygen exchange and some taxonomic responses varied. Sequence accession is ENA PRJEB8381.

## Eligibility and extraction inference

Neither component supplies independently replicated light units: field site is confounded with treatment, and incubation shelf with light. Neither12 lamps nor sediment cores become independent treatment replicates. Before/after measurements help describe changes but do not create spatial replication. The same cores pass through temperatures, so temperature levels are not independent experiments. Do not assign an inverse-variance field/light effect using the sample-level standard error.

This is the existing Westhavelland lighting infrastructure; do not count it as an independent lit landscape from insect studies there. Oxygen change can reflect net photosynthesis as well as respiration and may cross zero; do not force it into lnRR. Community composition and sequence relative abundance are not plant decline, pollination or measured predator attack.

The printed coordinate uses latitude minutes69, invalid sexagesimal notation. Do not silently reinterpret this as decimal latitude for remote-sensing moderators. A separately verified site coordinate is required. Retain the study as contextual experimental evidence with `hold_unreplicated_illumination_units`; no numerical effect added.
