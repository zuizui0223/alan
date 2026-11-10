# ALAN | 30 Korean islands, real moth trapping (2007–2025) × historically prior ALAN — actual empirical Results

2026-11-10, independently source-verified real field data. ALAN PR #35 (draft stacked on VIIRS500m PR #34).

## Source and chronology audit (completed, not inferred)

- Independently published Choi et al. (2026), Biodiversity Data Journal DOI 10.3897/BDJ.14.e211519, data DOI 10.5281/zenodo.21438943.
- All original published files individually Zenodo-MD5 verified, CP949 species taxonomy decoded. **30 real Korean coastal islands, 332 distinct original georeferenced light-trap locations, 377 distinct island × site × observed-year samples (17 actual years spanning 2007–2025), 1,179 identified moth species, 30,442 actual captured moth individuals**. No unknown years/empty samples or species labels invented. See original [sampling support ledger](../results/alan_real_korean_island_moth_sampling_support_20261110_verified.json), [source QA Actions](https://github.com/zuizui0223/alan/actions/runs/42772784877).
- Original 2025 source has **2 islands / 23 site-year observations only**, not the full 30 islands: assigning a 2024 500m raster to the older 2007–2024 moth samples as contemporaneous light would be anachronistic.
- For **all 377 genuine samples**, sampled exact coordinates from the **actual prior calendar year's Li & Cao corrected PCNL v3 annual original GeoTIFF** (survey minus one, actual raster years 2006–2024), verifying every originally published MD5. All 17 source satellite years were retrieved in 4 disjoint groups and original sample IDs preserved. **377/377 source site-year matches have valid PCNL** (0 missing); **58 original island×survey-year groups** have at least 2 true spatial sites with distinct satellite values before the more restrictive >=0.02 light-span gate. [Original-source 17-year join Actions SUCCESS](https://github.com/zuizui0223/alan/actions/runs/42773386737), artifact 13742604065; output SHA256 8a790f57024c93e486f957f15454339896ed9cfd967fff43ef77d24085a18a72.
- Original source does NOT provide precise month/night/duration or trap-night effort by site-year. Thus this analysis concerns **sampled light-trap capture assemblages**, not unbiased population density, individual moth breeding success, effective nocturnal pollen delivery or genuine plant reproduction. Light-trap phototaxis is an ascertainment mechanism.

## Main true ecological source fit

Prespecified exploratory **within actual island × observation year** WLS: response ~ log1p(previous-year corrected original PCNL pixel) + observed site elevation + island×year fixed effects. Inference clustered by actual independent islands, not across 30,442 individuals. Each actual island contributes equal total weight; genuine spatially supported groups only (>=2 sites and >=0.02 log1p prior light span). Every source-verified outcome is shown, not selected by nominal significance.

| Source-fitted observed light-trap endpoint | Genuinely supported islands / site-years | beta per log1p prior PCNL | 95% nominal island-cluster CI | nominal p |
| --- | --- | ---: | ---: | ---: |
| **Primary** expected richness among 20 moths sampled without replacement (finite rarefaction) | **27 / 239** | **−0.40194 species** | [−0.92778, +0.12390] | **.12822** |
| Rarefied 10 moths (support/sensitivity) | 27 / 292 | −0.20304 species | [−0.37168, −0.03439] | .02018 |
| log(1 + actual moths trapped) | 27 / 347 | −0.06628 | [−0.17444, +0.04188] | .21902 |
| log(1 + actually photographed? **NO**: actually trapped identified moth species) | 27 / 347 | −0.07234 | [−0.16088, +0.01620] | .10505 |
| Shannon exp(H) effective moth species from actual counts | 27 / 347 | −1.16423 | [−2.84278, +0.51433] | .16585 |

The original 20-individual primary has no resolved negative association at nominal 95% CI. There is a **consistent negative direction in sampled moth richness**; one sensitivity at 10 sampled individuals has an unadjusted interval excluding zero. The model does NOT support the simplistic positive prediction of larger numbers of light-trap moths near higher satellite nighttime-light pixels: measured trapped-individual coefficient is mildly negative and uncertain. But lamp attraction operates at a much finer spatial scale than PCNL and is not directly measured by this macro sampling design.

[Exact originally measured model result](../results/alan_real_korea_30island_moth_count_diversity_priorpcnl_verified_20261110.json) · [successful original 30-island ecological model Actions](https://github.com/zuizui0223/alan/actions/runs/42774729931).

## Genuine source robustness and multiplicity

All outcomes/sensitivities use exactly the same original year-verified source, without outcome-conditioned photo acquisition. [Source executed robustness Actions](https://github.com/zuizui0223/alan/actions/runs/42775157869), [complete signed robustness JSON](../results/alan_30island_true_moth_priorpcnl_leaveisland_FDR_eras_verified_20261110.json):

- For the five prelisted main ecological endpoints, **the 10-individual rarefied-richness nominal p=.020175 has Benjamini–Hochberg q=.10088 and Holm familywise p=.10088**. NO five-endpoint corrected result meets q or adjusted p < .05. The primary 20-individual p=.12822 remains unresolved.
- Leave-one-actual-island-out across the 30 source-island names (including 3 not informative to the primary within-year spatial analysis): **30/30 fitted directions remain negative** for 10-individual rarefied richness, 20-individual rarefied richness and trapped count. Rarefied10 coefficients range **−0.23563 to −0.15348**; in all 30 variants the nominal island-cluster 95% CI upper bounds remain negative. Primary rarefied20 coefficients range −0.48458 to −0.19775 but all source CIs include zero. The 3 islands not informative in the source-spatial model contribute no independent sign changes; do not advertise 30 independent statistical replications.
- Excluding earlier survey cohorts: 2012+ (27 supported islands, 252 rarefied10 samples) β=−0.19881 [−0.36898,−0.02864]; 2017+ (26 islands/225 samples) β=−0.25296 [−0.46764,−0.03829]; 2020+ (22 islands/141 samples) β=−0.36655 [−0.65079,−0.08230]. These are heavily overlapping post-hoc sensitivities and cannot substitute for a new independent external replication.
- For the 20-individual primary, site-year-equal vs equal-island weights changes p: island-equal β=−0.40194 [−0.92778,+0.12390], siteyear-equal β=−0.44920 [−0.88134,−0.01706] nominal p=.04220, islandyear-equal β=−0.46699 [−0.97291,+0.03894]. The direction persists but nominal statistical resolution **depends on weighting**.

## Independent biological interpretation and next controls

The results are an important step beyond GloBI plant-source labels: these are original moth specimens from real traps, under independently observed long-term satellite nighttime light. But **lighting can alter phototactic attraction to the trap, the plant/herbivore community, moth phenology and trap sampling**. Without recorded trap-night duration, survey dates, local weather and historical urban land use, one cannot separate these alternative mechanisms or infer absolute moth population decline.

This is a **30-island one-region moth data validation**, not replicated proof of a global nightlight effect, and not a proof that moonlight, moth-mediated pollination or flower-colour evolution were affected. Integrating this real moth source with the previous 2024 500m dark-corridor source requires a time-matched 500m annual VIIRS series (where supported) or an explicit **prior-year PCNL-derived corridor** with consistent resolution for each surveyed year. The two 2025 moth-survey islands alone are not an adequate global spatial effect model.

Recommended restrained research statement: **"Across 27 spatially supported Korean islands, higher historically prior satellite nighttime light was associated with lower rarefied moth richness in some exploratory specifications, but the primary 20-individual estimate and five-endpoint multiplicity correction do not establish a confirmatory effect; captured moth abundance had no resolved positive light association."**

Do NOT present the effect as demonstrated 'moths attracted to light increase local effective floral pollination', causal decline in moth density, plant-colour selection or global biodiversity loss.
