# Independent check of Ceppi2025 descriptive extraction

Checked directly from `raw_data_ceppi.xlsx`, sheet `Leaf_area_consumed`, on2026-11-20. This calculation did not import the main extraction's estimates or aggregate table. The original workbook SHA256 is `b7c0c7a81315db06951cabcfcb070f926236a42d1f05c0b93b65333cad88e78d`.

## Source-definition check

[Primary Methods2.1,2.4,2.5](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.70040) specify43 nights and21,35,70lux versus dark. Published severity models use only damaged leaves. The raw `4000_*` labels are retained as source labels; paper describes nominal3700K lamps. They are not independently calibrated colour-temperature measurements.

## Independent calculation

Restricted records to `control`, `4000_30`, `4000_50`, `4000_100`, leaving840 leaves:105 per site/treatment. The extra `2200_100` tree is excluded by treatment scope. `Perc_Herbivory` is complete. Damage means strictly positive `Perc_Herbivory`. For each installation, calculated (a) mean over all105 leaves, (b) damaged fraction, (c) mean over damaged leaves only. Then weighted the two sites equally. This is not a ratio of summed leaf areas or an equal-weight average of tree-specific conditional means.

| Endpoint | lux | RR | lnRR sampling variance |
|---|---:|---:|---:|
| All-leaf mean percentage |21|0.441825|0.048805785988|
| All-leaf mean percentage |35|1.149192|0.635822493501|
| All-leaf mean percentage |70|0.618031|0.001035522226|
| Damaged-leaf fraction |21|0.867470|0.064864657201|
| Damaged-leaf fraction |35|0.939759|0.012437382905|
| Damaged-leaf fraction |70|0.614458|0.013509325583|
| Mean percentage among damaged leaves |21|0.438354|0.133177974784|
| Mean percentage among damaged leaves |35|1.156937|0.744576286769|
| Mean percentage among damaged leaves |70|1.138948|0.002232308488|

For each contrast, the two independent-check scores were `lit_site / mean(lit_sites) - dark_site / mean(dark_sites)`. Covariance is their joint sample covariance divided by2. It has rank1; diagonal values are above. For example covariance(all-leaf21,all-leaf35)=0.176158498360, covariance(all-leaf21,all-leaf70)=-0.007109112191, and covariance(all-leaf70,conditional70)=0.001520396348.

The approximate t1 interval for all-leaf70lux lnRR is[-0.890096,-0.072337]; all other eight intervals include zero. The relatively small variance of this one contrast arises from cancellation in two matched-site scores. It does not supply robust precision across forests or establish a general intensity threshold. Nine contrasts remain dependent measurements from two sites; no new independent experiments are created. The conditional endpoint's different weighting and selected damaged leaves prohibit interpreting its difference from the unconditional endpoint as formal mediation.

After the extraction outputs became available, aligned by endpoint/treatment and checked all9 ratios,9 variances,18 interval bounds and81 covariance entries against this direct raw-data calculation: agreement at relative tolerance1e-12 (covariance absolute tolerance1e-15). The implementation's intervals are on the response-ratio scale; the interval above is on the log scale.
