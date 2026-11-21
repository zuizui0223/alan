# Penryn: four complete species, shared-unit covariance, trophic modification

2026-11-21. **One experiment, 54 physical mesocosms, four selected species;
not four independent experiments or a general wetness effect.**

The new analysis compares the ALAN biomass response of *Lotus pedunculatus*
with the geometric mean response of *Agrostis tenuis*, *Anthoxanthum odoratum*
and *Hypochaeris radicata*. All 216 species-by-mesocosm observations are
positive, and all 36 species-by-treatment survival counts agree with original
Supplement S1. *Centaurea nigra* has a complete biomass table but discrepant
survival counts and is excluded. Missing species observations are never zeros.

The primary author data are [Figshare 20304039](https://doi.org/10.6084/m9.figshare.20304039),
original [PDF](https://ndownloader.figshare.com/files/36265164), SHA256
`baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0`.
The original [Anic et al. 2022 paper](https://doi.org/10.1002/ece3.9157),
Methods 2.1, identifies six randomly allocated mesocosms per light-by-trophic
cell, each illuminated unit with its own LED strip. Original treatment codes
are L1 dark, L2 amber LED and L3 white LED; T1 plants only with repeated
insecticide/molluscicide, T2 plants plus herbivores and T3 plus predators.

## Estimand and covariance

Within each trophic condition, compute each species' lnRR from arithmetic
arm mean biomass. Subtract the mean of the three comparator lnRRs from the
Lotus lnRR. Exponentiation gives a **ratio of proportional biomass responses**,
not Lotus's share of total community biomass. The two light-specific marginal
results are equal-weight averages over the three trophic conditions on the
log scale; they are not ratios of pooled arithmetic means.

The 24 original species effects have a complete 24 by 24 sampling covariance.
For each participating treatment cell, centered normalized biomass values
contribute signed `(x / mean(x) - 1) / sqrt(n * (n - 1))` to a score column,
aligned by physical mesocosm. Its Gram matrix is the first-order delta
covariance. This retains both same-unit species covariance and shared dark
controls. Explicit linear weights give six cell contrasts, two marginal
contrasts and four trophic modifiers. Their full 12 by 12 covariance is saved.

Intervals are approximate t5 intervals, using six units per treatment cell
as the conservative degrees-of-freedom choice. Bonferroni intervals are
reported separately for two marginals, four modifiers and six cell contrasts.
This does not remove post hoc species selection or provide exact small-sample
coverage. Neither a p-value nor a meta-analysis diamond is added.

| Contrast | Ratio | Approximate nominal 95% interval | Family interval |
|---|---:|---:|---:|
| Relative response, amber | 0.635880 | 0.363211–1.113246 | 0.319206–1.266716 |
| Relative response, white | 0.491652 | 0.242467–0.996925 | 0.205994–1.173439 |
| T2 vs T1, amber | 0.936496 | 0.223109–3.930923 | 0.111722–7.850053 |
| T2 vs T1, white | 1.146607 | 0.210132–6.256572 | 0.092722–14.179080 |
| T3 vs T2, amber | 1.518673 | 0.498783–4.623985 | 0.291584–7.909777 |
| T3 vs T2, white | 1.010765 | 0.231005–4.422618 | 0.113382–9.010644 |

All two-colour family intervals and all trophic-modifier intervals include
one. The white nominal upper bound is only just below one; it is not a
confirmatory significant result after considering both colours.

## Trait support and sensitivity

The original Kew [BIFloraExplorer CSV](https://github.com/RBGKew/BIFloraExplorer/blob/master/data-raw/BI_MAIN.csv)
was independently downloaded and checked against Git blob
`e6a1da929af81f370d2a1b79177b489e3e573a06`. The four PLANTATT moisture F values
are 5, 6, 4 and 8 respectively; light L values are 6, 7, 8 and 7.
Agrostis tenuis is linked to A. capillaris through the
[Kew POWO synonym record](https://powo.science.kew.org/taxon/386415-1).
The printed biomass binomial has no authority, a limitation retained in the
trait receipt. These British/Irish habitat affinities are separate from the
European EIV source and from experimental soil moisture.

Only Lotus has F greater than 7. A moisture effect is inseparable from its
species identity and legume biology here. Comparators' average L is 7, as is
Lotus's L, but an ordinal habitat-light score does not control physiological
ALAN sensitivity. The same-L Lotus/Anthoxanthum comparison has ratios 0.750
and 0.508, with both nominal intervals including one.

All three single-comparator and all three leave-one-comparator-out marginal
point estimates retain the direction; this is not repeated statistical
significance. With one comparator at a time, ratios range 0.521–0.750
(amber) and 0.457–0.512 (white); all six nominal intervals include one.
A distinct mean-log-per-mesocosm estimand gives 0.665 and 0.515. Its point
estimates are saved separately, without claiming they are identical lnRRs.

## Files and reproduction

From repository root, with numpy, pandas and scipy:

```bash
python database/alan_experiments/plant_sources/analyze_complete_species.py \
  --out database/alan_experiments/plant_sources/complete_species_trophic
python -m unittest discover -s tests -p test_alan_penryn_complete_species.py
python database/alan_experiments/plant_sources/verify_complete_species_covariance.py \
  --analysis-dir database/alan_experiments/plant_sources/complete_species_trophic \
  --out /tmp/penryn_independent_verification.json
python database/alan_experiments/plant_sources/verify_complete_species_pdf.py \
  --pdf /path/to/anic2022_data.pdf \
  --selected database/alan_experiments/plant_sources/complete_species_trophic/selected_source_rows.csv \
  --out /tmp/penryn_original_pdf_check
python database/alan_experiments/plant_sources/verify_complete_species_traits.py \
  --source /path/to/Kew_BI_MAIN.csv \
  --out database/alan_experiments/plant_sources/complete_species_trophic
python database/alan_experiments/plant_sources/plot_complete_species.py \
  --source database/alan_experiments/plant_sources/complete_species_trophic \
  --out manuscript/biotic_synthesis/figures/penryn_complete_species_trophic
```

`species_selection.csv` records exclusions. `selected_source_rows.csv` retains
the selected original values and printed number strings. `species_effects.csv`
and `species_covariance.csv` give all 24 effects and 576 covariance entries.
`contrast_weights.csv`, `composition_and_trophic_contrasts.csv` and
`contrast_covariance.csv` specify every derived comparison.
`comparator_sensitivity.csv`, `mean_log_estimand_sensitivity.csv` and
`mesocosm_log_indices.csv` preserve the alternative comparisons.
`analysis_summary.json` fingerprints input files.

The analytical tables were generated with Python 3.12, numpy 2.3.5,
pandas 2.2.3 and scipy 1.17.0; the CI workflow pins those versions. SciPy
1.14.1 gives slightly different extreme Student-t inverse-CDF values: its
largest family-interval endpoint difference here is 1.8e-11, while effects
and covariance agree at the strict 1e-12 verifier tolerance. Using the
generation versions keeps that tolerance intact and all saved results unchanged.
Plot replay additionally requires matplotlib; PDF source verification requires
Poppler's `pdftotext`.

The independent source check used Poppler layout extraction and Decimal
comparison, not the original repository extraction function. It matched all
216 keys, 432 treatment labels, 648 biomass values and page references. The
independent numerical check used 36 arm means and a separate 24 by 36
Jacobian, matching all effects and covariances to less than 1e-14. Receipts
are saved here. These are AI checks, not human duplicate extraction.

## Scientific use

This is a targeted response comparison that complements the existing three
stage-specific exploratory meta-analyses. It does not increase the number of
independent programmes, identify predation or competition as a mediator,
estimate pollination or survival, or explain the macro recording-rate result.
T2 minus T1 is the complete treatment package, including the difference in
pesticide regime. T3 minus T2 evaluates predator addition; broad intervals
provide no precise evidence of its modifying effect.

Earlier 13-species PLANTATT associations on PR88 did not retain shared-unit
covariance and predate the Centaurea survival reconciliation. This analysis
answers a different, fixed-species comparison and does not replace those
historical diagnostics or the current 108-effect catalogue.

[Japanese results and figure](../../../../manuscript/biotic_synthesis/PENRYN_COMPLETE_SPECIES_ja.md).
