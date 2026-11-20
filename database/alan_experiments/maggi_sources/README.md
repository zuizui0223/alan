# Maggi2020: site-level original-data extraction

Paper: https://doi.org/10.1111/1365-2435.13485 . Data: https://doi.org/10.5061/dryad.6m905qfw2 ; recovered mirror: https://zenodo.org/records/5003041 . All three recovered workbooks match Dryad SHA256 digests. Source files remain in the external cache. See [independent workbook inspection](../search/maggi2020_data_check.md).

## Available endpoints and units

Four sites (two illuminated, two dark) are independent lighting units. Three quadrats per site and grazer treatment are averaged before estimating sampling variance. There is no pairing between dark and illuminated sites. Nine snail observations per quadrat are repeated measurements, not independent replicates. Final snail extraction uses time73. All exclusion-quadrat final snail counts are zero; no exclusion lnRR or pseudocount is used.

The pigment workbook contains carotenoid/chlorophyll ratio and scytonemin, not fluorescence biomass, quantum yield, or chlorophyll concentration. The bacterial workbook contains sequencing counts, not organism densities. Neither supplies the missing biomass response. Final pigment measurements and earlier fluorescence measurements are not interchangeable.

Within-file site identities are retained. Snail `light=l1` means dark; pigment `LightSite=l1` means an illuminated site. These labels must not be joined. Cross-workbook physical-site mapping and overlap with Maggi2018 remain unresolved, so covariance is calculated within the pigment workbook only.

## Calculation and result

Run `python database/alan_experiments/maggi_sources/extract_sites.py`. The script verifies original hashes, aggregates quadrats and produces site means, lighting effects, exclusion-minus-grazed contrasts, and covariance. Delta-method variance uses the two independent site means per arm and retains the covariance between grazer treatments and pigment endpoints within sites. An independent linearized-site calculation checks the contrast covariance. Approximate intervals use conservative t1 critical values; they are not exact small-sample confidence intervals.

The ratio of lighting response ratios (excluded/grazed) is0.6156 for carotenoid/chlorophyll (approximate95% interval0.3011–1.2584) and0.9062 for scytonemin (0.000641–1280.95). Neither excludes1. The final grazed-quadrat snail count response ratio is0.8214, with extreme uncertainty from two sites per arm; it is not a grazing or predation rate. Original full-time-series conclusions are a different estimand.

These are supplementary extraction products, not additions to the108-effect main catalogue. They do not test whether lighting removed a positive grazer effect on biomass: that endpoint is absent from the recovered data. Four pigment effects and two derived contrasts share the same sites; their numbers do not represent independent experiments. Site covariance has rank2 and must not be regularized merely to obtain a pooled model.
