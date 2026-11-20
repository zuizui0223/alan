# Deitsch2023: original-source audit and extraction hold

Original [published article](https://par.nsf.gov/servlets/purl/10490269), DOI10.1098/rspb.2023.0153; [publisher supplement](https://doi.org/10.6084/m9.figshare.22140093.v1). Cached source hashes and four-round calendar are in `database/alan_experiments/deitsch_sources/`.

Eight forest plots form four lit–dark pairs. LED exposure lasted7.5weeks. Plasticine caterpillars were deployed for72hours in four rounds; the fixed final-round endpoint is24–27July2021. Independent n=4 per arm, not521 recovered models or16 plot-rounds. Arthropod attack marks provide a standardized attack proxy. They do not measure actual adult-moth mortality, parasitism, or reduced plant reproduction. Bird and mammal marks require separate endpoints; uncertain parasitoid-like marks cannot establish parasitism.

The paper reports higher arthropod attack incidence under lighting. Its season-pooled summaries combine repeated rounds and cannot supply final-round paired variance. Predator and parasitoid catches vary with sampling method; activity-density is not attack rate.

## Access and what remains

The article identifies original data/code at [EDI](https://doi.org/10.6073/pasta/71525c9c724db0a21249584afca9ddbd), package`knb-lter-hbr.368.2`. Portal access returned Turnstile; no bypass attempted. Public PASTA metadata endpoint explicitly returned403 unauthorized. The separate publisher supplement was legally accessible and retrieved, but its nine tables contain model results, plot distances and taxon summaries rather than plot-round attack numerators and denominators. No original individual-data file was recovered.

**Zero effect sizes extracted.** No pooled-season SD was substituted, no missing attack was treated as zero, and no digitized approximation was presented as raw data. Need four paired final-round proportions, recovered-model denominators and identifiers. Then compute arm means/SDs and paired covariance across the four pairs; preserve all four rounds separately and estimate joint covariance only from matched units. Zero means require a non-lnRR endpoint decision rather than an arbitrary correction. Existing positive seasonal result remains narrative evidence pending valid extraction.
