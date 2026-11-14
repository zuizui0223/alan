# ALAN — multi-species nocturnal moth activity, corrected original Dryad source gate

**Original independent external evidence:** Briolat et al. (2026), Proceedings of the Royal Society B, [DOI 10.1098/rspb.2025.2704](https://doi.org/10.1098/rspb.2025.2704), source dataset [Dryad DOI 10.5061/dryad.jq2bvq8ps](https://datadryad.org/dataset/doi%3A10.5061/dryad.jq2bvq8ps), CC0 published metadata. This is an experimental animal behavioural activity mechanism, **not a direct conspecific stigma pollen, actual flower visitation or plant reproductive fitness study**, and should not be treated as 23 independent published ecological field studies.

## What is actually known, and in which source tier

The published paper reports:
- **843 independently sampled adult moths** from **23 species** across **three families** (Erebidae, Geometridae, Noctuidae) in its main experimental design; not 843 moth population lineages and not a global species census.
- **85% mean activity reduction** under white LEDs at **10 lux** relative to natural nighttime light. Other amber and RGB LEDs also suppress activity; lower illumination can still matter.
- Individual repeated video frames are not independent insect specimens, and activity does not establish actual flower visits or pollen delivery.

These are **published original authors' results, not newly refitted coefficients**.

## Genuine source metadata retrieval (executed)

[GitHub Actions #43570752664](https://github.com/zuizui0223/alan/actions/runs/43570752664) verified the public Dryad REST dataset and all corrected source file metadata. Its official latest original published dataset **versionNumber 4**, version ID 413461, was issued **2025-12-11** after an earlier 2025-12-03 version; the authors corrected a source movement-column issue, so the original R code must be followed. Public corrected files include:

| Original corrected file | Source compressed/actual file size | SHA256 from official Dryad manifest |
|---|---:|---|
| Activity_MainData.csv | **77,107,372 bytes** | e08798c04d5708e5ea51e0ceda8f4ea64f192b38656726c1f7b9a94305f18c04 |
| Analysis_Script.R | **79,730 bytes** | f33cd8af210527dbd98d8c278731f53fea0fefc01485ad0945f925478274586e |
| Species_size.csv | **764 bytes** | c7b872843166acbf55ca04df3662592f9bac7f5127250ea80730539b6d761e6e |

Main data corrected file ID **4507228**, analysis script file ID **4507224**, light spectra another source file. Public file lists and SHA are confirmed. Neither original frame rows nor the authors' exact source movement computation were retrieved in this stage.

**Direct source-byte fetch is blocked by Dryad authentication**: [Actions #43571082852](https://github.com/zuizui0223/alan/actions/runs/43571082852) returned **HTTP 401** when requesting the original R script bytes; it then emitted an honest *source metadata verified, raw access blocked* result (archive [artifact #14104159001](https://github.com/zuizui0223/alan/actions/runs/43571082852/artifacts/14104159001)). There is **no** computed per-moth activity reanalysis and **no** invented taxon-specific vulnerability effect. The script is set up to verify the exact real file bytes and count *moth IDs* rather than frames whenever authenticated source access is available. The workflow exits successfully on the observed evidence gate and does not imply that data bytes were downloaded.

## Why this fits the macro research programme

This source expands the **candidate mechanistic support** for ALAN responses of *nocturnal pollinator animal channels* from an individual species to three families and 23 species. But the measured outcome is activity suppression, not flower visiting, pollen transfer, compensatory fly/bat/gecko pollination, autonomous selfing, fruit/seed production or time-mismatch.

Compare with separate, genuine complementary evidence:
- Lu/Tibet (five blocks): real nocturnal Lepidoptera visits decrease, Diptera visits increase, pooled 11-plant moth-majority grouping **fails** reliable independent block transfer.
- Giavi/Switzerland (six experimental sites): daytime insect guild response signs are site and ecosystem dependent; even the Swiss pooled Diptera sign reverses from 5/6 individual site directions.
- Boom/Netherlands (four sites): same moth pollinates flowers and its offspring consume seeds; spectral lighting strongly alters infestation while no significant net seed-output response was resolved.
- Briolat/UK (23 moth species) supplies **actual nocturnal animal-activity suppression under controlled light**, but its source individual-grain results have not been independently reprocessed here due Dryad 401. These are four studies, with different response variables and nested designs: do not pool effect percentages as if they estimate the same functional syndrome×ALAN beta.

An ecological species-level moth-activity loss does **not** automatically imply a pollination-service loss. That mechanistic transition must be verified via pollen deposition and local fruit/seed outcomes, coupled with flowering phenology and the site's alternative animal visitors. Island vs mainland history and actual predation / floral service dual-role remain separate moderators to be measured rather than guessed.

## Reproducibility

Public original data-version file-discovery pipeline: scripts/analysis/alan_moth23_Dryad_original_source_manifest_20261114.py and test/workflow. Corrected data individual-grain source audit (fails closed on missing download rights): scripts/analysis/alan_moth23_Dryad_corrected_source_grain_20261114.py, matching test, original GitHub Actions workflows. The public results contain metadata hashes, HTTP status and source design only; no video frame data or fake reanalysis of 23 species. This work is a Draft stacked after #67 and leaves all previous ALAN findings and main unchanged.
