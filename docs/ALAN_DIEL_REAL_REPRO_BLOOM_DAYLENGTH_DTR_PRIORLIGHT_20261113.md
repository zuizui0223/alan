# ALAN: Does natural daylength and temperature range condition the functional value of night flowering?

**Original-source replay and environmental site/year matching: SUCCESS (2026-11-13).** [GitHub Actions #43418576694](https://github.com/zuizui0223/alan/actions/runs/43418576694), [aggregate JSON artifact #14031943793](https://github.com/zuizui0223/alan/actions/runs/43418576694/artifacts/14031943793), [Draft PR #56](https://github.com/zuizui0223/alan/pull/56).

## Original experimental source and limits

Original independent Kendall & Nicholson published field data, Ecology Letters 2025 DOI 10.1111/ele.70036. This original paper **already meta-regressed daylength and daily temperature range separately**. The additional analysis here is **joint field natural conditions × botanically documented bloom timing × historical (nonmanipulated) nightlight**, including strictly held-out botanical family/study/plant predictive transfer, **not a claim of discovering daylength or DTR associations for the first time**.

Parent source: 301 authentic original study×plant×field site×start year×month-season×fruit/seed outcome comparisons where night-only and day-only pollinator exclusions actually measure fruit/seed reproductive output. 118 independent publication studies, 123 plant species, actual field dates 1993–2022. Original all 30 Li–Cao corrected annual satellite PCNL 1992–2021 GeoTIFF MD5 verified at exact original published study field coordinates, mean of annual log1p strictly BEFORE field-study start year; actual exposure comparison is mainly between studies. Parent PR #55 original botanic blooming and direct flower-form phenotype complete-case population 196 source paired field outcomes, 73 independent studies, 74 species, 28 botanical families.

New original author independent environmental site-season sources: **218 distinct original reported study×field location×start year×start/end-month daylength entries**, Git blob 97ace3708bc356d87ce4fcd2e00fc4f1987be949, and **218 original daily temperature range entries**, Git blob f2f16164e96570e6acd5e6ab975295d2bcc6bd1f9c. Source daylength CSV uses old Latin-1 bytes; the pipeline verifies original Git SHA before decoding. DTR uses UTF-8. No fuzzy study matching: exact original normalized study ID, lat/lon each rounded 5 decimal places, original start year and original two month labels; conflicting or missing environment never imputed. Original source climatological/astronomical daylength is a calendar-season and latitude quantity, NOT remotely observed nightlight or field lux.

**Result:** of 196 original complete-floral-trait paired reproductive observations, 194 match exactly to original published daylength, 196 to DTR, **194 have BOTH** (72 independent studies, 73 plant species, 28 families). There are 48 source night-bloom studies and 28 day-bloom studies. Only 7 families contain both types in this aligned sample. Natural daylength ranges 10.66–18.99 hours, DTR 5.6–20.1 Celsius degrees. Original source daylength and DTR study-location/year/month key files have zero conflicting duplicates.

## 1. Same sample study-cluster regressions (reproductive FUNCTION, not visitation)

Response is the actual original source night/(night+day) fruit/seed outcome from genuine within-study pollinator-exclusion pairs. Study-equal weighted WLS and original study-cluster CIs, true study year, geographic polynomial, published elevation, fruit-versus-seed metric and authentic previous satellite PCNL control. An identical observation cohort is used whenever adding environment. Study×species two-way cluster uncertainty also calculated.

- Baseline bloom night vs day among the 194 matched pairs: **beta +0.19539 fraction** (+19.54 pp), 95%CI [+0.08415,+0.30664], p=.000576; study+species two-way p=.00261.
- Control for original published daylength+DTR together: night bloom **+0.25936 fraction** (+25.94 pp), CI [+0.15441,+0.36430], study-cluster p=.00000127, two-way p=.0000178. This is observational across source plants/studies.
- Add botanical FAMILY FE on the same 194 real pairs: bloom +0.21417 fraction (+21.42 pp), CI [+0.02725,+0.40108], p=.0247; study+species two-way p=.0284. The family FE only uses the 7 families with both bloom classes for identification, and may be fragile.
- Daylength main effect in the non-family-FE joint model +0.08449 fraction per day hour, p=.032 study-cluster but two-way p=.067; DTR main −0.01315 fraction per Celsius degree, p=.093. With family FE, DTR −0.03666 fraction per Celsius degree, CI [−.06080,−.01251], nominal p=.0029 and 2way p=.0105, **but** without family FE it is not supported. Do not call this robust causal thermal selection.
- Original **night-bloom × daylength** interaction alone is **−0.06480 fraction per natural daylight hour**, study-cluster CI [−.12090,−.00870], p=.0236, 2way p=.0358. Under simultaneous DTR interaction, night-bloom × daylength **−0.05869**, CI [−.11697,−.00042], p=.0484 study, **TWO-WAY p=.0646 with CI [−.1210,+.00366] crossing zero**. The interaction with DTR itself is not resolved (two-way p=.305).
- **Night-bloom × historical prior artificial nightlight after both natural-environment controls**: beta **−0.01092 fraction per mean log1p PCNL**, 95% study-cluster CI [−.07042,+.04857], p=.719, 2way p=.708. NO supported ALAN-driven loss of night-blooming reproductive advantage from these observational data.

All natural environmental terms may proxy different climatic zones, altitude, study taxa and geographic selection, rather than measured insect/pollen visit changes. Original climatological daylength is dependent on latitude and season; it is not an independent light-at-flower treatment.

## 2. True cross-family generalization on 194 identical source pairs

Frozen fixed-regularization Ridge alpha=5 with train-only feature scaling and per-TRAIN publication-equal fit weights; leave each of **28 botanical families** out completely and predict actual night reproductive share. Identical actual original 194 pairs, same geographic/year/elevation/fruit-seed/prior-PCNL controls, no original paper/species/family identity used as a feature. Each family has equal *test total* weighting in the LOFO-family table.

| Model | Held-out family RMSE (actual night reproductive share) | Improvement vs original field context |
|---|---:|---:|
| Original geography/year/PCNL/metric/elevation | **0.30213** | baseline |
| Add **night/day flowering** | **0.27385** | **+9.359%** |
| Add original **daylength and DTR only** | 0.30726 | **−1.699%** |
| Add flowering + DAYLENGTH only | **0.26165** | **+13.398%** |
| Add flowering + DTR only | 0.27109 | +10.273% |
| Add flowering + both daylength+DTR | 0.26426 | +12.535% |
| Add flowering × natural conditions interaction terms | 0.26501 | +12.284% |

**Key distinction:** daylength/DTR alone worsen out-of-family prediction relative to the original field context, but when combined with bloom time can improve RMSE beyond bloom alone. The largest numerical benefit in the family-holdout comparison comes from independently compiled seasonal/daylength information plus bloom timing, not adding an abundance of correlated flower-color/shape traits. Whether this is genuine ecological moderation rather than study selection/latent environment covariates requires independent validation.

### Original family-level uncertainty: can the extra climate data be said to improve over bloom ALONE?

The same 28 source botanical families were evaluated separately. Only **15/28** have a lower held-out family-level squared error from bloom+daylength+DTR versus bloom-only. Incremental family-equal RMSE improvement **+3.50%** versus bloom-only. A fixed-seed 3,000-replicate bootstrap of the original family blocks returns **DESCRIPTIVE 95% interval [−3.75%,+8.92%]**, spanning zero. **Not a confirmed universal climate-dependent night-pollination functional syndrome.**

## 3. Truly unseen published studies and species

5-fold block on entire published STUDY rather than random rows, with same 194 paired field results:

| Source model | Held-out study RMSE | Improvement vs original field context |
|---|---:|---:|
| Original context | 0.28703 | baseline |
| Bloom only | 0.27497 | +4.200% |
| Daylength+DTR without bloom | 0.28864 | −0.562% |
| Bloom+DTR | 0.26994 | +5.954% |
| Bloom+daylength | 0.26644 | +7.172% |
| **Bloom+daylength+DTR** | **0.26323** | **+8.291%** |
| Bloom+both environmental interaction terms | 0.26487 | +7.719% |

**Even stricter source:** Remove all plant species seen in each held-out published STUDY fold from its training source (12, 44, 21, 18, 21 train rows additionally purged per fold). Study/plant purged held-out RMSE: context **0.30152**, bloom-only **0.29421**, bloom+daylength+DTR **0.28213**, bloom+environment interactions 0.28355. Hence the bloom+environment model improves **~4.11% over bloom-only** when both entire studies and their plant species are unseen at training. This is a descriptive cross-validation sample, not independent environmental manipulation.

## 4. Scientific adjudication

**Supported as an exploratory macroecological *source pattern*:** flowering phenology strongly corresponds to actual measured diel reproductive allocation across independent botanical taxa, even after original daylength/DTR controls. Natural seasonal daylength and temperature range may contain conditional information that improves prediction of real original day/night fruit/seed function in a selected global study sample. The exact study/site/season climate matches 194/196 eligible originally measured reproductive cases, making the ecological question materially more credible than interpreting GloBI visitor labels or flower photographs as effective pollination.

**Not supported:** a causal photoperiod-induced shift in effective pollinators, an experimentally induced flower-timing effect, a universally transferable cross-family prediction gain, a confirmed independent natural-temperature mechanism, effective moth-specific pollen deposition, a nightlight-caused pollination change, compensatory autonomous selfing or a native island-only population syndrome. The original paper already tested environmental conditions and bloom traits; the new conditional pairing/transfer results are exploratory extensions, not a discovery of the base botanical effect. Low within-study opposite-bloom support and unrandomized source geography/timing remain.

**Next crucial mechanistic data:** paired real observed at-flower nocturnal/diurnal effective pollen deposition and independently quantified native community visitor availability on the same sites, with season/daylength/daily temperature range and actual local lamp lux; replicate across islands and across botanical lineages. Need independent study populations to distinguish environmental channel availability from taxonomy/observer choice.

## Code, restricted data, and audit chain

Source script: scripts/analysis/alan_original_daylength_DTR_diel_repro_function_macro_20261113.py; synthetic tests: tests/test_alan_original_daylength_DTR_diel_repro_function_macro_20261113.py; original SHA-source GitHub workflow .github/workflows/alan-real301-diel-daylength-DTR-original-year-20261113.yml. Original source / parent PR 55 unchanged; still Draft PR 56, main unchanged. Source inputs are original publicly published study environment data with both original 218-row season tables, source SHA original 301 authentic night/day reproductive source and all six annual satellite-PCNL shards from completed original source GitHub Actions. Public artifact contains **aggregate-only** results, no original study plant×field coordinate, rare taxon trait or site season rows.
