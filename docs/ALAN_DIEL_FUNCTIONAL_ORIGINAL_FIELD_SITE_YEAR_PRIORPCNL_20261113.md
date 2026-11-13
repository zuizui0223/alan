# ALAN: Global original matched daytime/nighttime reproductive experiments versus prior nightlight

**2026-11-13 independent author-source preflight SUCCESS**: [GitHub Actions #43379308366](https://github.com/zuizui0223/alan/actions/runs/43379308366). Research original: Kendall & Nicholson 2025, *Pollination across the diel cycle*, Ecology Letters DOI 10.1111/ele.70036. Original author Git study blob SHA1 39a82d966afcdde670603add1c881ace8791ffe9 and plant trait blob SHA1 18381e550a7e405fe4f198909f0bb381bbecf1e0. The published final 135 papers and 139 plant species are distinct from the 140 raw source study IDs.

## Scientific correction and true original source dimensions

Coarse original study × plant × metric source has 231 fruit/seed matched night/day treatment pairs, but several papers have multiple individual field sites, true experimental years and months. It is biologically incorrect to use one site for a response pooled across such sites.

The corrected source unit is original study × species × fruit/seed metric × original measurement unit × reported lat/lon × experiment start/end year × month bounds, with both actual daytime and nighttime treatment records and nonnegative numeric same-unit reproduction outcomes. The original raw treatment source of 1652 rows yields **301 valid originally site/year/method-matched reproductive comparisons** from **118 original independent studies, 123 published plant names, 135 real field coordinates and 163 study × site × experiment-start years**. Original field start years 1993 to 2022. These source metrics are **178 fruit-set** and **123 seed-set** outcomes; an author-study/plant/site may contribute both, so 301 is NOT the independent-study sample size.

270/301 paired outcomes involve plants described as native in their source; 197/301 coordinates were explicitly supplied in the original study text, 98 inferred from place names, 6 digitized from published maps. Original author elevation measurements are linked to all 301 locations, and only 171/301 have valid day/night treatment-duration information. The original experiment end year is not reported for **284** observations, with 17 explicitly multiyear start-to-end cases. The prior-light exposure uses only years before the **reported start**, never a guessed completion date or article publication year. Excluded 37 original site-period comparisons before 1993 (no earlier PCNL); 4 two-arm outcomes with both means zero; 1 invalid numeric treatment.

Original paired-source SHA256 21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af. Original author source files are retained in transient CI; source details need not be republished.

Descriptive pairing result, prior to reading any nightlight: median original night/(night+day) fractional fruit/seed response 0.583 and per-study-equal weighted mean 0.570. These are not the author manuscript's standardized-effects phylogenetic meta-analysis, not proof that night pollinators are generally better, and not actual moth-only pollination.

## Genuine original historical light linkage

Official original [Li and Cao corrected annual PCNL v3 Zenodo record 17013414](https://zenodo.org/records/17013414). The earlier original GloPL proof-of-source code from ALAN was copied unchanged to sample genuine annual calibrated 1992-2021 GeoTIFFs (30 years, six independent 5-year shards). Each raster has original Zenodo MD5 verification. Only PCNL from 1992 up to actual published field start year minus 1 contributes to each exposure mean; missing raster pixels are unknown, never coded unlit. This is an observational site light history, **not experimental manipulated lamp light**.

Source code:
- scripts/analysis/alan_original_DvN_fieldsite_year_source_20261113.py
- scripts/analysis/alan_authentic_diel_siteyear_PCNL_6shards_20261113.py
- scripts/analysis/alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113.py
- .github/workflows/alan-original-diel-field-priorPCNL30y-20261113.yml

[Original-source raster workflow](https://github.com/zuizui0223/alan/actions/workflows/alan-original-diel-field-priorPCNL30y-20261113.yml) fits study-equal study-clustered exploratory original night-vs-day experimental reproductive response against historical previous light after adjusting source study year, original author elevation, geography and fruit/seed measurement differences. Source-only quality restrictions predeclared independently of light effects: actual native status, coordinates provided in original article, original hours per period, original 1993-2022 study period.

**Never claim:** light caused a change in actual insect abundance or effective moth-only pollination, historical light exposure is an intervention, the published study species is geographically island-native without site verification, or 301 treatment comparisons are 301 independent species/populations. Previous ALAN PR #43-50 remain unchanged, PR #51 Draft on prior PR #50 and main unmodified.
