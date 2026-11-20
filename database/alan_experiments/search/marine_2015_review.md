# Davies2015 marine settlement: primary Methods and Results reviewed

Reviewed 2026-11-20. [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4424621/?pdf=1), DOI [10.1098/rsbl.2015.0080](https://doi.org/10.1098/rsbl.2015.0080).

## Source-derived design and findings

At Menai Strait, 36 panels occupied 18 randomized boards for 12 weeks from July 2013. Six boards received each treatment: control, 19 lux or 30 lux. Two panels per board were summed for counts and averaged for cover. Forty-seven taxa occurred; 13 met the authors' frequency threshold. Cover composition differed among treatments; count composition did not unless one outlying replicate was removed. Species responses included increases and decreases. The Methods state lights switched on at dawn and off at dusk, inconsistent with the stated nocturnal treatment; retain this unresolved wording.

Supporting data are listed as `rsbl20150080supp3.csv`. The PMC attachment request returned HTML rather than CSV; original measurements are not yet acquired.

## Extraction inference

Use boards, not panels or taxa, as replication. Both illuminated arms share six controls. Preserve the outlier in any raw-data primary extraction and report author outlier sensitivity separately. Do not select only the five significant taxa. Colonization/counts/cover are separate endpoints; none directly measures predation, pollination or plant reproduction. Nighttime operation needs corroboration before final exposure adjudication. Hold effect extraction, not full-text reading.

The [Bolton2017 primary paper](https://doi.org/10.1016/j.scitotenv.2016.10.037), Methods 2.1 and 2.3, places its experiment in Sydney Harbour with 2014 plates. This excludes identical sites/plates with Menai2013; it is not a complete overlap audit of every coastal paper.

## Access record

Requests retrieved the full PMC HTML (140,225 bytes) via the `?pdf=1` URL, despite a browser challenge at the bare URL. Local cache: `../gbif_annual/alan_experiment_sources/marine2015/article.html` and `article.txt`. EuropePMC XML timed out. Supplement endpoint: `https://pmc.ncbi.nlm.nih.gov/articles/instance/4424621/bin/rsbl20150080supp3.csv`; its 200 response was HTML (1,817 bytes), not data. Do not treat HTTP success as successful CSV retrieval.
