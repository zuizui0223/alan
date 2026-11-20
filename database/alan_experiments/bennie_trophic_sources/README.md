# Penryn trophic follow-up: slug abundance

Paper: https://doi.org/10.1111/1365-2664.13240 . Original data: https://doi.org/10.5061/dryad.8217n76 . Mirror: https://zenodo.org/records/4988926 . All13 files match original Dryad digests. Source workbooks remain in the external cache.

## Shared experiment and endpoint

This is the Penryn54 programme already used by Bennie2015 and Anic2022, not a new experiment. All36 slug mesocosms match existing treatment allocation. All54 rows of the2013 final aphid endpoint match the previous extraction, including missing entries; these are not new independent observations. Slug counts are not attack rates. The predator treatment combines ground beetles and ladybirds, rather than isolating a single predator species.

The paper restricts slug inference to observations through September2014 because of a population collapse and2015 reintroduction. We keep the final September2014 endpoint: both dark-arm means are zero, so all four final lnRRs are undefined. Counts and empirical mean-difference variances are saved, without pseudocounts. An observed zero sample variance is not proof of zero population uncertainty and is not an infinite-weight input for meta-analysis.

## Separate secondary estimand

Source versions are not identical throughout2013: across378 matched aphid plot-dates, B6 on2013-08-29 is7 in the recovered2018 table versus0 in the2015 source, and C4 is4 versus0 on that date. `aphid_2013_source_discrepancies.csv` preserves both values. Final2013 values agree; no prior result is overwritten and no intermediate source is silently preferred.

The supplementary period analysis averages15 supplied monthly counts within each mesocosm: June–November2013 and January–September2014. December2013 is absent and is not zero-filled. Equal weighting is by recorded month, not a time-integrated population measure. Each treatment cell has six independent mesocosms, not90 replicates. We do not pair by plot letter.

Four lighting lnRRs and two predator-presence contrasts use arm means of these mesocosm averages. Shared dark-arm covariance is retained. Conservative delta-method t5 intervals are approximate. The predator-contrast ratio is1.111 (0.512–2.412) for amber and0.672 (0.232–1.945) for white. Neither interval excludes1. White lighting's mean-count ratio is0.905 without predators and0.609 with predators. This period summary differs from the author's longitudinal GLMM; it does not refute or reproduce that model's reported test.

Run `python database/alan_experiments/bennie_trophic_sources/extract_slugs.py`. `audit.json` records identities and exclusions; `programme_crosswalk.csv` anchors physical units; final-zero and period results remain separate. These effects stay outside the108-row main catalogue. Any later joint analysis with Penryn aphids/plants must reconstruct their shared-mesocosm covariance rather than count the papers as separate experiments.
