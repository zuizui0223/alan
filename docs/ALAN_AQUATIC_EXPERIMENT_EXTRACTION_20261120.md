# Aquatic predator primary-source audit, 20 November 2026

## Conclusions

No new controlled spider-abundance effect is ready for the main meta-analysis. Meyer2013 is an uncontrolled before/after manipulation; Sullivan2019 has a promising separate wetland experiment but raw download failed; Ogbeide2026 has verified independent replication and downloaded spider isotope measurements, not abundance.

## Meyer & Sullivan 2013: important negative evidence with a different denominator

[Primary author PDF](https://bpb-us-w2.wpmucdn.com/u.osu.edu/dist/1/24733/files/2020/03/12-2007.1.pdf), DOI10.1890/12-2007.1; [original supplements](https://doi.org/10.6084/m9.figshare.c.3295991).

Methods pp1323–1324: observational sampling covered nine reaches (three per ambient-light class), June2010–June2011. In July2011 five reaches were illuminated continuously at10–12lux: three previously low and two previously moderate (one moderate reach lost access). ALL five experimental reaches were illuminated; there was explicitly no unmanipulated2011 control. August2010 versus August2011 paired comparison therefore confounds added light with year. Spider webs within1m of stream and up to2m height were an abundance surrogate. Three emergence/pan traps per reach are subsamples, not independent light replicates.

The reported44% spider decline is valid as the authors' before/after result, but cannot be entered as contemporaneous lit/dark lnRR. AppendixD mean/SD values describe ambient-light groups by month, not experimental2011 lit/control. For example August observational Tetragnathidae(no./reach): high31.9 SD18.1, moderate45.1 SD14.3, low36.4 SD10.2. No effect was calculated. This protects the contrary result without misclassifying design. Year-adjustment cannot be recovered from an absent control.

## Sullivan, Hossler & Meyer 2019

[Paper](https://doi.org/10.1002/eap.1821), [original Dryad deposition](https://doi.org/10.5061/dryad.gb9534p).

Dryad explicitly separates UrbanStreamCategories_2010, UrbanStreamGradient_2010-2011, and WetlandFieldExperiment_2013. Therefore the entire paper must NOT be classified purely observational; the wetland component remains eligible for design review. Stream years/authors/geography overlap Meyer2013 and must be crosswalked before independent study counting. Dataset filenames alone do not establish exact individual-record duplication.

The six primary files and their metadata were identified (version8787). Public API metadata succeeded; API downloads returned401 and website file_stream downloads403. Error responses have been renamed .download_error.json in external cache; they are NOT data. Wetland experiment raw file38946 and README38947 are highest-priority access targets. No n, means, SD, or allocation were inferred from abstract, no new effect extracted. The abstract's negative orb-web abundance claim cannot yet be assigned specifically to the experimental component.

## Ogbeide et al2026: independent aquatic interaction experiment

[Paper](https://doi.org/10.1111/1365-2435.70335), [original data/code](https://doi.org/10.5281/zenodo.14516568).

Primary methods sections2.3–2.6 verify16 artificial streams, randomized fully crossed ALAN×signal-crayfish, four independent mesocosms/treatment,24May–5July2022. Three emergence traps per flume and three female Tetragnatha extensa selected for isotope analysis per flume are subsamples. Spider consumer isotope sampling occurred week6. Thus48 spider isotope records mean n16 units, NOT48 independent replicates.

Downloaded mixture_spider.csv and source_spider.csv and associated primary isotope files. Derived flume means first, then group mean and sampleSD among four flumes. These summaries are in aquatic_sources/ogbeide2026_spider_isotope_*.csv. They are an audit of the raw-source replication, NOT new main meta-effects: isotope ratios are signed baseline-dependent quantities; lnRR is inappropriate. ALAN-alone/control d15N means8.056667/7.932500, among-flumeSD0.116078/0.118239; d13C−26.041667/−26.008333, SD0.332103/0.172530. Aquatic diet proportions require the authors' mixing model and baseline uncertainty; these isotope differences do not establish reduced aquatic prey or spider abundance.

The322MB code archive is identified but not downloaded. Companion2026 DOI10.1186/s12302-026-01345-4 appears to concern this same experimental programme; exact units/dates must be verified before counting independent evidence. No diet credible intervals were converted into samplingSD.

## Implication for wetland mechanism

Evidence is compatible with changes in prey routing and taxon-specific predator responses, not a uniform wetland increase in spiders. None of these sources yet measures a chain from adult-moth predation to pollen limitation to wetland-plant decline. Wetness is directly aquatic/riparian for these experimental settings, but coordinate-derived distance/water-frequency covariates remain uncomputed. Maintain outcome-specific strata and independent-study grouping; do not pool before/after, isotope, and predator-abundance denominators.
