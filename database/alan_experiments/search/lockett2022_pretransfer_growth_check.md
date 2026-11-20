# Lockett2022 pre-transfer growth check

Independent raw-data calculation, 2026-11-20. This is an explicitly supplemental pre-transfer endpoint, not terminal biomass or evidence of a trophic mechanism.

## Primary-source summary

[Methods](https://doi.org/10.1016/j.envpol.2022.119803) describe diameter measurements before lighting at week 1 and again at weeks 11, 18 and 38. Psyllid introductions and plot transfers began later. Measurements were taken near 300 mm height, with an exception for short seedlings, using two perpendicular caliper readings. Trees in poor condition remained physically in plots but were excluded from the authors' experimental analyses. Poor condition includes discoloured foliage, stunting or visible invertebrate impacts; it is not synonymous with death. Lighting followed a checkerboard arrangement, not randomized allocation. Twelve plots per arm, each initially containing six trees, are the lighting units. This design can support a pre-transfer plot-level contrast, with systematic allocation and possible spatial dependence retained as limitations.

## Independent original-data checks

Source: [Mendeley original tree workbook and README](https://data.mendeley.com/datasets/pfrrkgs878/1), `Growth` sheet. There are 144 distinct trees at each of four weeks. All 288 week-1/week-18 diameters are positive and observed. Each initial plot contains six trees with one consistent lighting label. No biomass-sheet treatment labels were used. The crossed-plot problem found at harvest does not apply before transfer.

`diam_change` equals `diameter / diamweek1`, not the difference or the ratio minus one: maximum absolute discrepancy across nonmissing observations is 4.97e-9. Week-1 change entries are missing. The new calculation uses the two original diameter columns, not these derived changes. README calls diameters centimetres, whereas Methods state caliper resolution in millimetres without providing a definitive archived-unit conversion. Keep source values as recorded; a common unit conversion cancels in the growth-ratio estimand.

### Condition flags and interpretation

The author flag `use4weeks1-18` excludes six trees: light tree 41, and dark trees 78, 87, 100, 121 and 124. All six have a poor-condition flag at week 11. Tree 121 is flagged at weeks 1, 11 and 18; tree 124 at 1 and 11; the other four only at 11. At week 18 their recorded diameters are respectively 5.90, 3.45, 7.52, 6.75, 1.30 and 7.52. Positive stem diameter does not establish living tissue, health or growth in biomass.

These records do not preclude a descriptive **full-cohort recorded stem-diameter change** analysis: no values are fabricated, missing measurements imputed, or poor trees silently excluded. They do preclude calling it living-population growth, survival-adjusted biomass, or randomized intention-to-treat evidence. Preserve flags and identify selection differences. The author-retained cohort is shown below only as an inclusion-rule diagnostic, not to select the more favourable result.

## Independently calculated effect

For each original plot, calculate the arithmetic mean diameter separately at weeks 1 and 18. Within each arm, average the 12 plot means equally. Estimate

`log(mean_D18_light / mean_D1_light) - log(mean_D18_dark / mean_D1_dark)`.

For each arm use the empirical 2×2 covariance of plot means divided by 12 and gradient `(-1/mean_D1, 1/mean_D18)`. Sum arm variances; retain within-plot time covariance. No spatial pairing is inferred. The t11 interval is a conservative working interval, not protection against arbitrary spatial correlation or systematic allocation bias.

| Cohort | Light trees | Dark trees | ln contrast | Variance | Ratio | t11 95% interval |
|---|---:|---:|---:|---:|---:|---|
| All recorded trees | 72 | 72 | −0.035041265561 | 0.007482724362 | 0.965565570831 | 0.79817092–1.16806670 |
| Author growth flag only, diagnostic | 71 | 67 | −0.038031703737 | 0.008081081551 | 0.962682419777 | 0.78986842–1.17330611 |

Full-cohort arm means (week 1, week 18) are light (1.4586111111, 5.6365277778) and dark (1.5051388889, 6.02375). Their log growth ratios are 1.351783544493 and 1.386824810054; variances 0.004375322184 and 0.003107402178. The author-flag diagnostic retains all 24 plots and equal plot weighting despite unequal numbers of retained trees.

Verdict: suitable as a transparent supplemental descriptive plot-level extraction under the stated working independence assumption. It provides no clear directional treatment contrast and does not resolve harvest label errors, identify plant survival effects, or demonstrate herbivory mediation. Do not promote it to the prespecified final endpoint or choose it because terminal extraction is currently blocked.
