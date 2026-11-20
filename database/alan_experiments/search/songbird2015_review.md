# de Jong2015: field reproduction, shared lighting programme

Reviewed 2026-11-20. Full [primary Methods and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC4375368/), DOI [10.1098/rstb.2014.0128](https://doi.org/10.1098/rstb.2014.0128).

Seven Dutch sites contributed in2013, eight in2014. Each site had four randomized lighting transects (white, green, red, dark), five posts per transect and nine nest-boxes per transect. Two songbird species were studied. Light effects on great-tit laying date and chick mass occurred in2013 but not2014; most breeding and return outcomes showed no significant lighting effect. Models included site, with nest-box additionally used for chick mass. Fledgling counts excluded complete brood failures; great-tit brood failure was analysed separately. Adult survival was scored from return to breeding between years. The setup is the Dutch Light on Nature infrastructure.

The [institutional dataset record](https://pure.knaw.nl/portal/nl/datasets/data-from-effects-of-nocturnal-illumination-on-life-history-decis/) points to [Marine Data Archive file VLIZ_00000444_5821cc7422ac7](https://mda.vliz.be/directlink.php?fid=VLIZ_00000444_5821cc7422ac7). The archive returned a request/login page, not data. No request was submitted.

## Extraction inference

Recommend **eligible_field_design_pending_raw_data_and_programme_crosswalk**. Lighting replication is transect within site, not lamp, nest-box, chick or brood. Repeated sites across years do not create fifteen independent sites; shared controls, species, broods and years require covariance or explicitly selected endpoints. Individual settlement is not randomized, despite randomized illumination transects.

Final fledgling counts conditional on successful broods cannot stand alone for overall reproductive output; failures must remain represented. Return-to-breed measures apparent survival and can include dispersal/detection. Lay dates require day differences, not lnRR of arbitrary calendar values. No measured moth predation, prey abundance, pollen or plant response is supplied here.

Crosswalk sites/transects against other Light on Nature papers before pooling. Programme overlap is established, but exact overlap with each later study is not yet checked. Do not substitute model F-tests or figure SEs based on broods for independent-transect effect variances. No effect extracted. The request-page response is cached externally at `alan_experiment_sources/songbird2015/archive_response.bin`.
