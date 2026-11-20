# Parkinson 2023: original capture data and shoreline placement contrast

Status: original trap counts extracted; biomass held; no master catalogue changes. Extraction date 2026-11-20. Not independently re-extracted by a second reviewer.

## Primary evidence and scope

[Article DOI 10.1098/rstb.2022.0361](https://pmc.ncbi.nlm.nih.gov/articles/PMC10613537/), Methods 2(a–c), Figure 3; [original Figshare collection](https://doi.org/10.6084/m9.figshare.c.6856475). Methods establish 30 randomized lamp/pan-trap stations: three per colour × habitat combination. Treatments were green, violet, red, white, and dark. The experiment ran one night, 3 June 2022, at Waterfowl Pond, Michigan. Littoral traps floated over water; riparian traps rested on adjacent soil. Stations were approximately 5 m apart; habitats approximately 15 m apart. Thus the habitat contrast is shoreline placement, not wetland versus dry upland. Aquatic/terrestrial taxon origin, coded in Figure 3, is a separate variable.

Illuminance averaged 40.8, 13.7, 43.1, 42.5, and 0.4 lux respectively; dark readings included dusk light. Violet was not illuminance matched. Captures describe local attraction/interception, not population decline, emergence production, predation rates, or pollination. Independent treatment assignment does not replicate ponds or dates or eliminate possible shared-source attraction effects.

## Reconciliation of original spreadsheet

Source file `rstb20220361_si_004.xlsx`, SHA256 `1e39be30482475d04643cce942a568ee6e99814ac2773deeb43ab792b03aef86`. Retrieved through Figshare article API 24212174; raw files remain outside the repository in `../gbif_annual/alan_experiment_sources/parkinson2023/`. Original download metadata records URL, file ID and MD5; downloaded files matched those MD5 values.

Directly summing the `raw site data` sheet yields **4,699 captures**, agreeing with the published total. Thirty blocks begin where habitat and colour are populated; a stray `22` in a separator's date cell does not begin a block. The five dominant families total **4,475** in raw records (95.23% of captures). However, `overall abundance ` sums to **4,465**, and `top5families` counts sum to **4,448**. These are not interchangeable whole-community endpoints.

`cross_sheet_discrepancies.csv` records the exact conflicts. At raw Excel row 47, the five-family sum is 120 but summary abundance is 110. Four further discrepancies concern Cicadellidae entries in dark-control rows of `top5families`: raw values 27, 0, 5, 2 versus sheet values 0, 5, 2, 0. This pattern suggests a spreadsheet alignment issue but does not establish its cause. The raw sheet is retained without correction. Biomass sheets are held pending verification of their family-count inputs and calibration mapping. Published NMDS IDs are saved only as unverified cross-sheet identifiers; primary keys are raw-sheet row numbers.

## Reproducible extraction

Run `python database/alan_experiments/parkinson2023_sources/extract.py` from the repository root. Outputs include 30 trap units, taxon counts, 40 contrasts, 784 covariance entries for 28 computable lnRR contrasts, and 12 computable habitat contrasts with their 144 covariance entries. All cells contain three stations. Endpoints: all captures, aquatic-origin captures, terrestrial-origin captures, spider captures, and Cicadellidae captures. The latter is a named family, not an exhaustive herbivore guild. Source spellings and the article's origin assignments are preserved, including potentially coarse family-level classifications.

No pseudocounts: all spider controls and littoral Cicadellidae controls have mean zero, so 12 corresponding lnRR contrasts remain missing. These raw observations can support a future count model, but are not finite log-response ratios. Spider presence in traps does not measure their consumption of moths.

For finite endpoints, lnRR = log(mean lit / mean dark), with ordinary delta-method sampling variance. Covariance includes shared dark stations across colours and observed cross-endpoint covariance within the same treatment stations. Littoral and riparian stations are unpaired. The covariance matrix passes diagonal/variance equality and positive-semidefinite checks (minimum eigenvalue approximately −6 × 10⁻¹⁶, numerical rounding). Small n = 3 makes uncertainty estimates exploratory.

## Within-study findings

For all captures, littoral minus riparian lnRR is +0.728 (green), +2.375 (violet), −0.234 (red), and +0.939 (white). These four contrasts share controls and must not be treated as independent. Aquatic-origin capture differences are +0.064, +1.732, −0.285, and +0.483 respectively. Therefore the clearest shoreline amplification is under violet, while red shows the opposite point-estimate direction. Spectral and irradiance differences prevent attributing this solely to wetness.

This supplies direct experimental evidence about redistribution of potential aquatic–terrestrial resources, with substantial spectral dependence. It does **not** demonstrate that wetland plants decline through spider predation, moth pollination loss, or a generally greater ALAN effect in wet versus dry ecosystems. Both sampled habitats belong to the same freshwater shoreline system; retain that distinction in any meta-regression.
