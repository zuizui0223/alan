# ALAN: Real night/day reproductive function × GloBI guild evidence — source identification audit

**Real original-source replay SUCCESS, 13 November 2026.** [Final GitHub Actions #43447139278](https://github.com/zuizui0223/alan/actions/runs/43447139278), [aggregate and DOI-only queue artifact #14049265688](https://github.com/zuizui0223/alan/actions/runs/43447139278/artifacts/14049265688). Draft PR #57. Prior branches and main unchanged.

## Source and response

Original Kendall & Nicholson 2025 *Ecology Letters*, DOI 10.1111/ele.70036: **301 genuine study×plant×field site×true year×season×fruit/seed reproductive-outcome pairs** under separate daytime-only versus nighttime-only pollination treatments, from **118 published studies and 123 original study plant taxa**. Original source cohort and 30 Li–Cao corrected prior satellite PCNL raster years (1992–2021) reused from SHA/MD5-verified workflows, never future-of-experiment radiance or actual light manipulation. The measured endpoint is fruit/seed reproductive night share; NOT actual effective moth pollination.

Original complete GloBI v0.11 interpreted stable source: **24,453,954** interaction rows, original original-source SHA256 **670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3**. Original published flower-interaction scanner identified **64,316** butterfly/non-butterfly Lepidoptera positive flower-use claims and **4,547 unique original positive plant names**. These are names, not independently verified effective pollination or site-level censuses. Original source scan [#42924993614](https://github.com/zuizui0223/alan/actions/runs/42924993614).

## True source overlap: only plant names, no study SITE × YEAR animal records

Among the 123 original field-reproductive plant names, 23 had an exact GloBI original plant-name match to any Lepidoptera positive flower use, contributing 63 of the 301 original measured day/night reproductive study-site/year/metric pairs.

| Source-positive animal guild flower-use category | Exact plant names | Actual experimental reproductive pairs | Original field published studies |
|---|---:|---:|---:|
| Moth-family positive at any global source location | 14 | 52 | 22 |
| Butterfly-family positive at any global source location | 15 | 39 | 22 |
| Both GloBI moth and butterfly positive | 6 | 28 | 13 |
| Moth positive, butterfly NOT DOCUMENTED | 8 | 24 | 9 |
| Butterfly positive, moth NOT DOCUMENTED | 9 | 11 | 9 |
| Moth positive, source excludes iNaturalist and trophiCH | 11 | 36 | 16 |
| Butterfly positive, same exclusions | 8 | 22 | 11 |
| Strict non-iNaturalist flower-visitation relation moth | 3 | 11 | 6 |
| Strict non-iNaturalist flower-visitation relation butterfly | 4 | 13 | 7 |

Full-source GloBI v0.11 moth/butterfly qualified flower interactions have **ZERO original source rows with BOTH georeference AND event observation date available**. Hence no GloBI record can be assigned to the SAME observed original plant fruit/seed experiment SITE AND YEAR. There is no evidence of local diurnal pollinator compensation or effective nocturnal moth pollen delivery from this source. Record absence means database omission, NOT absence of a biologically active pollinator.

### False apparent source-guild functional differences

Original field night-share MEDIANS measured per distinct published plant species, but grouped by independent source documented GloBI claims:

| Source filter | Moth positive / butterfly not documented | Butterfly positive / moth not documented |
|---|---|---|
| All original GloBI v0.11 | 8 plants, median night share 0.6138 | 9 plants, 0.5006 |
| Excluding iNaturalist and trophiCH | 7 plants, 0.4976 | 4 plants, 0.6611 |
| Strict observational flower-visit relation | 2 plants, 0.5351 | 3 plants, 0.5036 |

The direction visibly **REVERSES** after source-provenance restriction. These are DESCRIPTIVE source documentation medians; no inferential p-value is reported. A legitimate positive-guild class comparison was pre-gated on 20 distinct source plants and 25 independent field published studies PER positive class; all source comparisons fail this gate. Non-butterfly Lepidoptera is not synonymous with genuinely nocturnal adult moth behavior, nor do butterfly records measure all daytime plant visitors.

### Original artifact digest correction

The historical ALAN note quoted compressed positive-flag GZIP SHA256 ending ...a67d513d9bad74ece8e01, whereas the actual original GloBI full-scan Actions artifact has SHA256 **cf750192b4d641a2d904f4f74d5b87d41fc97e04822a67c513d9bad74ece8e01**. The archived ORIGINAL compressed file and receipt were independently retrieved and validated (4,547 names, original full 24.45m release checksum, 64,316 qualified flower claims). The initial workflow correctly refused the mismatched pin; corrected to physical archived original file SHA and successfully reran. Decompressed CSV flags independent digest: **20701811ca06d51f4a3e5f604c6e85cafa8673bd18875ed26591ff79738f6c5b**. Neither the original GloBI source nor previous research outcomes were silently rewritten.

## More promising: the authors' OWN field-pollinator observations

The 1,652-row original author data explicitly includes the field Pollinator_observations_included. Across the 140 original source study IDs: 111 papers mark some observations included, 27 do not, 2 source-unknown or mixed.

Among **118** published original papers supplying the **301** genuine paired fruit/seed outcomes, **92 papers indicate that some field pollinator observations were included**, representing **247/301 actual reproductive outcome pairs and 110 original plant species**. Of these 92 papers, **91 have their actual original bibliographic DOI**. This is a MUCH larger field-study-level acquisition opportunity than attempting to infer visitor identity from globally recorded GloBI plant names. However, the study-level observation flag does NOT specify animal identity, time, actual site or pollen delivery.

Author observation-positive field studies with plants that also have any global GloBI moth-positive relationship: 13 papers, 33 real day/night fruit-seed experimental pairs, 12 original source plant names. With independent stricter non-iNaturalist/trophiCH moth-positive plant claims: only 8 papers. GloBI positive global species and same paper observation inclusion do not establish local guild use.

### New public bibliographic metadata queue — not a confirmed guild matrix

Original-source Actions artifact additionally contains **ALAN_original_92study_public_DOI_pollinator_observation_curation_queue_20261113.json**: 92 original study identifiers, 91 original study DOI literals, original paper's number of reproducible field day/night reproductive outcome records, explicitly labeled unverified keyword screens of author MECHANISM notes and a curation status that is always **NOT verified local guild**. No original plant species identities, exact study coordinates, site-year plant records or original comments are included.

The queue is actionable for evidence curation. Each paper must be checked for actually observed moth, butterfly, bee, hoverfly, bat, bird, gecko etc. visitor guild at the SAME plant population, field site, true year and time window; original visit observation effort and effective pollen transfer must be recorded separately. Do not convert narrative suggestions or absence of a reported guild into biological observations. Papers with only study-level visitor observations but no matching trial/site/year must remain paper-level metadata.

## Overall adjudication

**Verified:** actual original reproductive experiments can be linked to global species-level flower-visitor documented guild claims, but only 23/123 original plants have any GloBI moth/butterfly flower evidence and 0 source claims have joint event georeference plus date. Source guild median signals reverse with provenance and fail independent study/species support. Original **92 published papers with 247 experimentally measured reproduction outcomes and author-reported pollinator observation inclusion** offer a more realistic basis for a new original site/time-guild functional database, and 91 DOIs are now queued as bibliographic metadata.

**Not established:** local nocturnal moth dependence, true pollen-transfer efficacy, light-induced pollinator loss, compensatory daytime services, island-only global ecological decline or adaptive pollinator shifts.

**Next reproducible data gate:** retain independent paper/plant/site/date/visitor ID and measured observation effort; first quantify actual site-year visitor datasets from the 92 source publications, require at least 20 botanically independent plants and 25 original published source studies in each genuinely confirmed and supported guild comparison, separately track source geographic and site annual ALAN. Continue source-neutral hypothesis tests regardless of sign rather than chasing nominal p-values from GloBI omission.

## Code / source provenance

- Original day/night fruit-seed × original GloBI source overlap script: scripts/analysis/alan_genuine_diel_GloBI_v011_flower_guild_overlap_20261113.py
- Original publication observation coverage and 92-DOI public metadata screening queue script: scripts/analysis/alan_original_diel_study_visitor_observation_coverage_20261113.py
- Synthetic unit tests: tests/test_alan_genuine_diel_GloBI_v011_flower_guild_overlap_20261113.py
- True authenticated GitHub Actions: .github/workflows/alan-original-diel-GloBI-v011-positive-guild-support-20261113.yml
- Final SUCCESS [#43447139278](https://github.com/zuizui0223/alan/actions/runs/43447139278) and aggregate+DOI queue [artifact #14049265688](https://github.com/zuizui0223/alan/actions/runs/43447139278/artifacts/14049265688). Draft PR #57, no preceding PR or main changes.
