# Zone-level monoculture biomass intervals

2026-11-20. This extends the existing four effects with uncertainty intervals;
it adds no study, experiment or catalogue effect. Reproduce with
`python database/alan_experiments/xue_sources/compute_monoculture_intervals.py`.

The original ScienceDB workbook was independently regrouped: 77 cylinders,
three ALAN and three dark zones, 2–4 surviving cylinders per zone/species.
All four zone-weighted means, lnRR estimates and delta variances agreed with
the existing extraction. Independent calculation here means a second analyst
agent, not two human extractors. The existing 4-by-4 covariance remains valid
and must accompany any future joint synthesis.

| Species | Light/dark biomass ratio | Approximate 95% interval | Welch df |
|---|---:|---:|---:|
| Hydrilla verticillata |1.0573|0.7402–1.5103|3.8295|
| Myriophyllum aquaticum |1.5719|0.7104–3.4784|2.3454|
| Myriophyllum spicatum |0.8393|0.4596–1.5326|2.5556|
| Vallisneria natans |1.0839|0.8909–1.3188|3.5924|

For each arm, v = sample variance/(n × mean²). The log-ratio variance is
v_light + v_dark; Welch df is (v_light + v_dark)² divided by
[v_light²/(n_light−1) + v_dark²/(n_dark−1)]. The interval uses a two-sided
t quantile before exponentiation. These are small-sample delta approximations,
not exact intervals; no multiplicity adjustment or pooled species test is made.
Each includes one. This does not establish equivalence or absence of effects.

The response is surviving-cylinder dry biomass under monoculture. It is not
survival, population decline, snail-mediated herbivory or an isolated causal
direct effect. The community experiment still lacks verified bucket-to-lamp
mapping and remains excluded from effect-variance synthesis. Signed feeding
mass differences remain unsuitable for an unmodified consumption lnRR.

Source: [article](https://doi.org/10.1002/ecy.70458),
[original data](https://doi.org/10.57760/sciencedb.35715).
The current local `xue2026/paper.txt` is a Cloudflare challenge response, not
a preserved full article; it must not be cited as a locally verified Methods
copy. The historical publisher reading is documented in the extraction note.
The preserved preprint is separately identified in `version_source_manifest.csv`.
This audit does not replace the final article with its preprint or certify
completion of a second human source extraction.

Independent follow-up read the preserved preprint Methods, lines127–137 and
182–190: herbivores were excluded; each cylinder held three seedlings of one
species; 80 cylinders were installed and 77 retained after mortality; zone was
an author-model random effect. This verifies the conditional scope in that
version, not a fresh read of final published Methods. The author square-root
response model is not the same estimand as the untransformed mean lnRR here.
All four Welch intervals were independently recomputed from the original
workbook and matched the displayed values.
