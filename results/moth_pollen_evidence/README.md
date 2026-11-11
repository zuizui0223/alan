# Moth pollen evidence among 544 shared plants

The frozen five-year comparison contains 10,698 site-years from 78 five-degree regions. Species labels derived from field and experimental literature apply anywhere a species occurs in this analysis. They do not assert that a moth carried pollen at the photo site or under the observed light exposure.

## Evidence

The species ledger accounts for all 544 accepted names. Positive evidence requires a species-resolved observation of pollen on a moth, species-level pollen DNA from sampled moths, or experimental evidence of pollen delivery through an identified moth visit. A plant occurring only in a study's diurnal columns does not qualify. Unresolved genus or aggregate records were excluded. Literal frozen names and unambiguous WCVP synonym resolutions were accepted; the source spelling Iris psuedacorus was corrected to Iris pseudacorus. The Spanish supplement reports certain plants at subspecies rank; these are represented by their parent species in the analysis. Percentages describe the screened evidence, not exhaustive coverage of each species' literature. Many taxa are represented by a single UK study. One survey of urban moths and bees (Ellis et al., 10.1111/ele.14261) could not be added from its linked data archive; combined bee/moth plant lists were not assigned to moths.

* Walton et al., [supplement S3](https://doi.org/10.6084/m9.figshare.c.4955030): positive counts in the seven nocturnal moth family columns.
* Macgregor et al., [Oxfordshire F4 pollen table](https://doi.org/10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708) and [Yorkshire pollen metabarcoding](https://doi.org/10.5281/zenodo.1408559): named positive pollen types in F4; species assignments with at least 50 reads in a named moth sample per the study method.
* [Spanish mountain networks, supplementary tables S1–S3](https://doi.org/10.1038/s41598-023-49944-y): positive nocturnal interactions inferred from three or more pollen grains.
* [New Zealand moth network incidence](https://datastore.landcareresearch.co.nz/dataset/moth-pollen-load-and-experimental-exclusion-data) and [Florida Table 2](https://doi.org/10.1111/een.13446) screened; no additional species matched.
* Species-specific evidence is referenced by source and type in pollen_evidence.csv. Direct stigmatic grain counts and single-visit reproductive outcomes retain separate labels. In particular the 43 positive species cannot all be classified as effectively moth-pollinated.

## Global light comparison

The same denominator (all photographed species in each site-year), row set, land balancing, covariates, region-year effects and cluster uncertainty were used throughout. The difference in slopes for the externally documented group was −0.172 pp per log1p radiance, and for the other 501 shared species it was −1.835. The two add to −2.007; adding +0.106 for the plants not shared between land classes reproduces the frozen −1.902. Their individual confidence intervals are in summary.json. The sizes of additive composition contributions are not tests of differential light sensitivity per species.

## Pollen limitation and self-compatibility

[GloPL authors' case table](https://github.com/idiv-biodiversity/pollen-limitation-data-descriptor/blob/abbe193770981b2c61713e75a0f01e004bdded8e/Output/GloPL.csv) contains natural-versus-hand-pollinated reproductive comparisons. The plant ledger records the unweighted median published log response ratio across available cases for each matched plant. A positive median indicates the observed direction and does not constitute a species-level significance test. The case table retains the reproductive variable, treatment, year and study citation. The historical experiments were not performed in the five frozen years.

[Self-compatibility compilation](https://doi.org/10.1093/aob/mcae056), Supplementary Data, Regression data sheet: use the source categories for SC, SI, self-sterile (SS), and partial compatibility; retain conflicts between multiple entries under an accepted name. The source includes some inferred SI mechanisms. No inference from family membership or lack of data was added. SC records show capacity to reproduce with self pollen; determining autonomous self-pollination or the fraction of reproduction supplied by moths requires separate observations.

Run scripts/analysis/alan_moth_pollen_evidence.py with the v0.11 input directory, signed MOIS directory, source-file directory and result directory. The study inputs are identified by source above and exact file hashes in summary.json. The script checks the frozen sample and numerical additivity against the existing model results.
