# Middendorp2025 ant foraging review

[Primary article](https://doi.org/10.1016/j.biocon.2025.111439), Methods 2.1–2.4, Results 3.2 and supplementary Tables S1–S2. Source URLs and hashes are in `source_manifest.csv`.

Twenty plots in a Western Australian woodland received alternating lit/control assignments (ten per arm). Four garden lamps illuminated each treatment plot. The experiment lasted 15 days; daytime and nighttime surveys were repeated. Five bait stations per plot offered honey and dead insect tissue. The authors reported no significant lighting effect on community or foraging responses. Plot-level data are available on request, not supplied as an open table.

The supplement contains five response models. Table S2 labels 28 August and 12 September as day, whereas the main Methods call them night; the other two dates are reversed. Table S1 also has ambiguous reference-level wording. Its protein-removal response uses a beta-logit model, whereas abundance, richness and incidence are log-transformed Gaussian responses.

## Extraction decision

We preserve 15 reported coefficients without converting them into arm-level lnRR. No arithmetic treatment/control means, independent-plot dispersion or joint coefficient covariance have been recovered. The group and interaction coefficients must not be added with an invented zero covariance. The plus/minus entries are transcribed as reported; no independently verified sampling variance is assigned.

Repeated photographs, bait stations and dates do not increase illumination replication. The four date/period conflicts are saved verbatim. These are reasons to hold quantitative inclusion, not grounds to discard the reported null result. Resource removal from dead bait does not measure live-insect predation, moth mortality or pollination.

Run `python database/alan_experiments/middendorp_sources/extract_reported.py`. The 15 entries remain outside the 66-effect catalogue.
