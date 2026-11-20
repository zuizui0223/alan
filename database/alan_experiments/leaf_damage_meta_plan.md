# Exploratory leaf-damage synthesis: specification before pooled calculation

Scope: retained, original-data-derived terminal leaf damage in Guo2025
(Hong Kong), Cieraad2023 (Netherlands), and Ceppi2025 (Switzerland).
This is a supplementary synthesis while systematic screening and programme
adjudication remain incomplete, not the final wet-versus-dry meta-analysis.

Include all available unconditional leaf-area-damage contrasts in these three
datasets. Exclude spider counts, incidence, severity conditional on damaged
leaves, chemistry and growth from this endpoint. Selection is not based on
direction or significance. Preserve plants lost before retained observations
as a limitation; no missing-to-zero replacement.

1. Cieraad: all six existing species-by-colour contrasts and full covariance.
   Two forest blocks; posts/trees are subsamples. No other Light on Nature
   report enters this synthesis.
2. Ceppi: all three all-leaf mean-percent contrasts and their joint covariance.
   Two sites. The other spectral treatment remains out of scope as previously
   documented; it has only one additional tree at one site. No other Swiss
   ALANeX report enters this synthesis.
3. Guo: both final-round plant endpoints. Average paired plots within each of
   four sites and treatment first; weight sites equally. Compute joint
   paired-site delta covariance. This differs from the prior seven-pair
   descriptive estimator and retains the broadest recorded clustering.
   Keep Aporosa at genus level because source species names conflict.

Within each programme, average log ratios equally over its fixed available
species/lighting menu, and use a'Va for sampling variance. This estimates
the average log response of that menu, not a response standardized to a
common species, dose or spectrum. Programmes are not weighted by how many
contrasts they contain. The three independent programmes then enter a
REML random-effects synthesis with modified Knapp-Hartung uncertainty
(variance inflation max(1,Q/(k-1)), t with k-1 degrees of freedom).
Retain all results regardless of interval width/sign.

Record source hashes, component weights, programme estimates and pooled
estimates. Verify with independent calculation. Small within-study cluster
counts (2,2,4) mean delta variances are uncertain even with a small-k meta
interval. Do not call tau2=0 proof of homogeneity, infer dose thresholds,
or infer predation/chemical mediation from damage. No moisture moderator,
species wetness effect, pollination decline or macro-plant mediation is
estimated. This synthesis does not alter the wet-affinity primary result.
