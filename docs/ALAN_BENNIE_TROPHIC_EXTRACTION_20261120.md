# Bennie 2015 original trophic experiment — extraction

2026-11-20. [Primary article](https://doi.org/10.1098/rstb.2014.0131); [official original supplements](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4375371/supplementaryFiles). Downloaded workbook `rstb20140131supp1.xlsx`, SHA-256 `009ea860f7493cb3ccd8eb0fa12503702f30eed1c38b2b570c4e0e6eb695f412`; accompanying methods/model PDF also inspected. Existing meta-analysis values were not used.

## Verified design and overlap

Fifty-four independently illuminated outdoor mesocosms crossed three light regimes (dark, amber, white) with three trophic communities, six mesocosms per cell. T1 contains plants, T2 plants plus introduced aphids/slugs, T3 additionally introduced ladybirds/ground beetles. Aphids are measured by standardized three-minute counts; *Lotus pedunculatus* flower heads by counts. Flower counts are a resource/phenology endpoint, **not floral biomass, fruit set or effective pollination**. Introduced predators constitute a treatment; their attacks were not quantified.

All **54 original plot IDs and both treatment factors exactly match** the existing Anic2022 biomass dataset. Verified mappings are Dark=L1, Amber=L2, White=L3; T1/T2/T3 agree for every plot. Therefore both reports belong to **one programme, `Penryn2012_2016`**: these 2013 counts and the 2016 final biomass harvest are different times/endpoints in the same mesocosms. They must not be counted as independent studies. The crosswalk is saved.

## Final-observation rule

All repeated source data are preserved: 918 plot-date rows, including 126 blank aphid observations in T1. Those blanks remain missing/structurally unmeasured, not fabricated zero counts. The final original dates are 12 September 2013 for aphids and 29 August 2013 for flowers. All relevant final treatment cells have six observed mesocosms.

| Aphids, 12 September | Lit mean | Dark mean | lnRR | Delta variance |
|---|---:|---:|---:|---:|
| Amber, T2 without predators | 21.3333 | 24.0000 | −0.117783 | 0.493329 |
| White, T2 without predators | 44.3333 | 24.0000 | +0.613683 | 0.473899 |
| Amber, T3 with predators | 36.5000 | 19.3333 | +0.635482 | 0.822453 |
| White, T3 with predators | 13.6667 | 19.3333 | −0.346871 | 0.993136 |

These four original-scale effects are extractable, with source means, SD and independent n stored, pending independent extraction verification. No transformation of counts or continuity correction was used.

For final flower heads, every dark-cell mean is zero. T2/T3 are zero under all lights; T1 has amber mean 0.1667 and white mean 0.6667. **All six flower lnRRs are undefined**. We retain original means and mean differences, flag zero means, and do not substitute an earlier date to obtain a favorable ratio. This final-season floor cannot represent the full flowering-season resource effect.

## Predator modification of the final aphid response

The contrast is ALAN lnRR with predators minus ALAN lnRR without predators, using four independent treatment cells per contrast:

- Amber: +0.753265; variance 1.315782.
- White: −0.960554; variance 1.467034.

Both are imprecise; directions differ between spectra. They do not establish that ALAN consistently increases predator control. The article's repeated-season model comparison likewise did not favor adding predators to the most parsimonious model. That published model selection result and these final-date contrasts are different estimands and should be reported separately.

The article's evidence for a resource-mediated plant–aphid pathway uses seasonal observations. The extracted final snapshot cannot on its own establish mediation from fewer flowers to fewer aphids. In particular, final aphid abundance does not show uniform reduction and final flower abundance is largely at zero. If seasonal flower supply becomes a target estimand, define a season-integrated count/phenology measure explicitly before analysis rather than silently replacing the fixed final-time rule.

## Dependence handled explicitly

`programme_covariance_including_anic2022.csv` contains the joint delta covariance for the four finite 2013 aphid effects and six existing 2016 biomass effects. Within an endpoint, spectra share the same dark mesocosms. Across years/endpoints, covariance is estimated from the six **same plot IDs in each shared treatment cell**, including shared control and shared lit cells; it is not assumed zero. Cross-trophic cells use separate randomly allocated mesocosms. With six units per cell these covariance estimates are uncertain and are not a substitute for programme-level dependence in a later meta-model.

The two predator-interaction estimates also share controls across spectra; their 2×2 covariance is exported. Undefined flower lnRRs never enter a covariance matrix or inverse-variance calculation.

## Files and checks

Own outputs are under `database/alan_experiments/bennie_sources/`: reproducible extraction, full repeated observations, final original units, arm summaries, effect estimates, predator contrasts/covariance, programme covariance, Anic crosswalk and hashed source audit. Checks verified 54 unique units per source sheet, six units per relevant treatment cell, complete relevant aphid data, unchanged source hash, and the 54/54 cross-year match.

Master tables were not changed; no commit. Four finite effect estimates represent **zero additional independent programmes** beyond Anic2022, not four new studies. The wet-affinity macro result remains unchanged.
