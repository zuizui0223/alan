# Maggi2020: recovered workbook coverage check

2026-11-20. Independently inspected every sheet with openpyxl, both formulas and cached values, and checked comments, hidden sheets and defined names. All three local SHA256 hashes match the [Dryad file metadata](https://datadryad.org/api/v2/versions/44480/files).

| Workbook | Actual data | Coverage |
|---|---|---|
| data_Littorinid_snails.xlsx |216 rows;24 quadrats×9 times; final time73| Snail counts, no direct grazing rate |
| data_Pigments.xlsx |24 rows;four sites×two herbivore levels×three quadrats| CarChl_ratio and Scytonemin only |
| data_Bacteria_Family.xlsx |32 sample rows in Foglio1; Foglio2/Foglio3 empty| Taxonomic sequence counts and read totals; includes two CA samples/site |

No inspected sheet, comment, named range or hidden sheet contains minimum fluorescence Fo, maximum fluorescence Fm, quantum yield or a separate biomass measurement. The pigment workbook contains a carotenoid/chlorophyll ratio, not the numerator and denominator as separate concentrations. **CarChl_ratio cannot be relabelled as photosynthetic biomass, chlorophyll concentration or quantum yield.** Bacterial read counts are not biomass either.

## Site and sample identifiers

- Snails: `light=l1` means dark, sites`s1,s2`; `light=l2` means lit, sites`s3,s4`. The explicit `lightherb` labels establish this coding. Each site has three control and three exclusion quadrats; no procedural-control counts appear. `plot=c1` etc. recur across sites, so plot alone is not a unique ID.
- Pigments: `Light=dark`, `Site=1,2`, `LightSite=d1,d2`; `Light=lit`, `Site=1,2`, `LightSite=l1,l2`. Thus the token `l1` has different meanings across files and is not a join key. Site number alone also repeats between lighting arms.
- Bacteria: `light=dark/light` crossed with `site=1,2`; plot prefixes C1/C2/L1/L2, with Control/Exclusion/CA and within-treatment replicate numbers.

Two dark and two lit sites are confirmed within each file. A cross-file map such as snail`s3`→pigment`l1` is plausible but not explicitly documented in the inspected files. Do not assume cross-endpoint covariance, or align individual quadrats across files by replicate number, without resolving this mapping. Within-file site aggregation is possible; quadrats and dates do not increase independent lighting replication above2/arm.

## Primary measurement distinction

[Methods](https://doi.org/10.1111/1365-2435.13485) distinguish in-situ fluorescence on18May from rock-chip pigment/bacterial sampling on26May. Fluorescence was used as a biomass proxy, separately from pigment ratios. Therefore this archive supports separate snail and pigment endpoints, while the reported biomass interaction requires another original source or documented figure extraction. The recovered files do not resolve the earlier experiment's shared-site status.
