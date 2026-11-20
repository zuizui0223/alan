# Cronin: frog signal receivers and predatory bats — version-specific audit

Reviewed 2026-11-20. This is primary-source extraction from **Cronin 2024 thesis chapter 4**, not verification of the final 2025 journal tables. No finite group-mean effect was added to the catalogue or existing pooled estimates.

## Sources and access

- Final article: https://doi.org/10.1093/jeb/voaf092 . Publisher web retrieval failed with a redirect error; the Groningen-hosted final PDF returned 403. Abstract and indexed excerpts are available, but final Methods/Results have not been fully read.
- Thesis: https://doi.org/10.5463/thesis.758 ; public PDF https://research.vu.nl/ws/files/352665541/acroninthesis%20-%2066cdf88acef78.pdf . Download succeeded (166 pages). Chapter 4 is labelled **In preparation**, PDF p.66. Methods/Results read at printed pp.69–77, PDF pp.70–78. Same title, four authors, taxa, experiment year and location establish a related report of the same experimental programme, not an additional independent experiment. Numerical identity with the final article remains unverified.
- Dryad: https://doi.org/10.5061/dryad.stqjq2cf5 . Public metadata and version 376568 file inventory recovered. Direct API downloads for code.R, data.csv and README.md returned HTTP401; no payload was retained as data and no alternative protected route was attempted. Web-rendered README was read. Reported repository digests are not locally verified file digests.

Raw accessible sources remain outside Git at `../gbif_annual/alan_experiment_sources/cronin2025/`. The manifest records downloaded-source hashes. Thesis PDF/text/images are not redistributed.

## Experiment and unit

Panama, August–November 2019: six locations separated by at least 750m, four sites per location (minimum site separation 25m), two paired call-choice plots per site. A single flashlight illuminated both choice plots at approximately 4.27 lux. Control, light, noise, and combined treatments rotated over four visits, separated by at least seven days. Monitoring was approximately four hours/night. Site, location and visit are repeated/grouped units; two speakers, approaching animals and video frames are not independent lighting replicates.

Responses are female frog approaches, female call choice, and fringe-lipped bat (Trachops cirrhosus) approaches to synthetic Engystomops pustulosus calls. Bat observations were classified using proximity and orientation/hovering/landing criteria, then presence/absence was modelled because count models had convergence problems. The bat model retained location but dropped site and visit random effects. These are not successful prey captures, mortality, fertilization, seed production or moth predation.

## Extracted original coefficients

`cronin_thesis_coefficients.csv` preserves the Light and Light×Noise rows from the interaction models in Table 2 (printed p.75) and Table 3 (p.76). All six are dependent coefficients from one experiment. Light coefficients apply to the reference noise condition, not an average across both noise levels. Do not pool log odds ratios with log count-response ratios. The thesis's prose describes a log link for choice while Table 3 explicitly says logit; the extraction records the table scale provisionally.

The female approach Light coefficient is −0.637 (SE0.220), bat presence −3.107 (SE1.117), and female choice −1.292 (SE0.466). Bat Light×Noise is −17.946 with SE724.078: this is an unstable fitted coefficient, not a precise synergistic effect. No back-transformed interval is constructed for it. Covariances needed to derive the effect of light when noise is present were not obtained. The interaction-free models are alternative analyses of the same data, not additional effects for pooling.

Further reconciliation needed: four visits ×24 sites imply 96 site-visits whereas Fig.3 states 99 observations; Dryad README lists visit numbers 1–5. Methods describe 293 females before exclusions and 282 approaches thereafter; counting and exclusion units need raw-data verification. The stated preferences and some back-transformations in prose do not consistently match exponentiation of table coefficients. Do not silently repair the text or infer missing raw cells.

## Decision and relevance

Keep final article as **partial primary/access hold**, not methods-and-results-read. Track thesis separately as grey literature, with overlap flagged. The three directions are relevant counterevidence to a universal increase of predators under ALAN, but they do not identify a wetland-plant decline mechanism. Do not combine this call-response assay with spider web prey stocks or aphid counts as a common predation-rate effect. No wet/dry habitat comparison, water proximity calculation or island contrast is supplied here.

A second automated reviewer independently confirmed the three Light coefficients, table scales, design nesting and major inconsistencies from the thesis text. This is not a second human extraction or a raw-data reproduction. Table 2 was also visually checked against the PDF rendering.
