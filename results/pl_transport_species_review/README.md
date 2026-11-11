# Species-level review of moth pollen transport

The review on 2026-11-11 covered all 53 plants initially classed as unconfirmed among the 64 shared plants with GloPL estimates. Searches used each accepted binomial and terms for moths, hawkmoths, pollen loads, pollen transfer and pollinia. The ledger lists one source and the reason for each decision. The original 11 positive classifications were carried forward. Primary observations were eligible whether made in the field or under controlled visits. A report of species-identifiable pollen on a moth, an attached pollinarium or measured single-visit deposition qualified. The Moorea field report is included as original data. Its publication and sampling context differs from the journal studies.

Unconfirmed means the material located in this search did not establish the link. Some searches yielded only abstracts or incomplete supplements. In particular the Portuguese nocturnal transport matrix (Banza et al., [Dryad](https://doi.org/10.5061/dryad.6ts275n)) could not be verified. Manning and Cutler's lowbush blueberry field record is a near match: the authors found tetrads on moths but documented other simultaneously flowering Ericaceae with similar tetrads. Neither blueberry nor the neighbouring rhododendron is promoted on those grains.

The model uses the original 78 regions and 10,698 site-years, including site-years without a group member. Species have equal weight in the regression of species-specific island-minus-mainland light coefficients on median pollen limitation. The model contrasts retain covariance within spatial clusters. A separate calculation fitting the species responses with explicit region-year indicators and regressing their light interactions on pollen limitation reproduced both updated slopes to within 1e-12. The pollen evidence can come from outside the focal study region and years.

From the repository root, supply the same two input directories as in the original analysis:

```
python scripts/analysis/alan_mating_light_contrasts.py INPUTS SIGNED results/pl_transport_species_review --review results/pl_transport_species_review/review_53.csv
```
