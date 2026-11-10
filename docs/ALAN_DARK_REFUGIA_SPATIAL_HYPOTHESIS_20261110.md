# ALAN | Spatial dark refugia and local connectivity: first real-source satellite pilot

**Status (2026-11-10): spatial exposure method implemented, with synthetic geometry tests. Real source metrics must be read from the GitHub Actions output; do not label synthetic tests as observations or invent measured island counts before execution.**

## Ecological question

Does the **spatial configuration** of dark land within an island influence the observed response to nightlight, after controlling the amount of dark land and the local mean light signal?

The strong biological hypothesis is that landscape-connected dark refugia might preserve movement and effective nocturnal pollination under artificial light. That mechanism **is not tested** by this first, purely geospatial source audit. No independent nocturnal pollination, pollen deposition, fruiting, fecundity or light-at-flower measurements enter this stage.

This is a deliberate change from a disputed post-hoc white floral ROI regression and from the original GloBI-documented moth-visit fraction. In the current source-audited 93m photo sample, the latter is profoundly sensitive to plant genus, popularity and the moth/butterfly study-claim record. Source-label differences are not a moth population effect.

## Source and eligibility frozen before evaluating any new response

- Input: exactly 85,318 original accepted-plant island 0.05-degree photographed cell×year records from the 93,209,351-photo provenance audit, joined to 1992–photo-year−1 source-verified PCNL. Validate the existing signed joined CSV against its original receipt and exact earlier-source SHA. Select actual 2025 photo-site records with complete historical source.
- Initial raster: **Li & Cao PCNL v3 2024** original 30-arcsecond (~1-km) intercalibrated DMSP-OLS/VIIRS yearly image. It **precedes** the photo year, not same-year VIIRS and not the 1992–2024 historical mean. Original DOI https://doi.org/10.5281/zenodo.17013414. MD5 enforced against the project's frozen official source manifest.
- Geography: original high-resolution **GSHHG 2.3.7** island unit and water polygons with original source ZIP/extracted SHA checked. Land-only mask from L1 island **minus L2 lakes**. Never treat sea pixels as dark habitats.
- Fixed **5 km radius geodesic** neighbourhood around each actual original photo-site coordinate, or a pilot random-less response-blind SHA-sorted maximum 8 sites per island in the 2025 year. The latter is a **capped diagnostic pilot**, not a probability sample of a global botanical population. Site IDs and observed year are preserved. (The hash order is computational sampling only, not ecological assignment.)
- Before looking at plant data, select thresholds PCNL native value <= **0, 1 and 3**. These are sensitivity settings in intercalibrated PCNL scale, **NOT calibrated lumens, lux or nW/cm²/sr at insect eye level**. No data-driven optimisation or significance-based threshold choice.
- Raster no-data, out-of-coverage, antimeridian sites and <16 original island-valid pixels, or <80% real raster land coverage, yield explicit unavailable status. Connectivity of zero dark pixels is **undefined**, not zero. Circle and pixels use actual source WGS84 geodesic distances; approximate pixel areas include cos(latitude). Four-neighbour edges define one connected dark component.

## Environmental variables actually measured

- Dark land fraction: dark original *land* pixels divided by all valid original island land pixels in the actual buffer.
- Largest dark component / dark land: a patch-continuity index conditional on any dark pixels; it is **not independent** of dark fraction and must be compared at similar dark coverage.
- Largest dark component / all land, four-connected dark patch count and distance from photo site to nearest dark raster pixel (right-censored at radius if none).
- Local original mean 2024 PCNL brightness, original land pixel count, source raster-coverage rate.

The project script outputs per **unique island × cell × 2025 × threshold** measures plus a signed source receipt and failure counts. Repeated thresholds on the same photographed site are NOT independent replicates.

## The later biological test, only when actual independent outcome support is demonstrated

Prespecify photo-site level within-island-year comparisons:

    observed_biological_outcome ~ log1p(local_PCNL) + dark_fraction
        + f(dark_fraction) + conditional_dark_connectivity
        + urbanization + landcover + climate + photo/observer effort
        + island fixed effects + observation-year fixed effects

Compare connected vs fragmented dark land **within overlapping bins of dark_fraction and mean local PCNL**. A sparse island with zero variation in dark connectivity is not a support island. Assess 5/10 km and each frozen threshold before interpreting a mechanism. Group cluster at the **island**, spatial-block sensitivity, leave-one-island-out, threshold sensitivity and explicit power/positivity denominators. If required, use a land-mask-constrained spatial permutation null with the same dark pixel count and buffer geometry; conditional index alone does not prove movement corridors.

**Do not fit this as an effect of pollination on current GloBI positive data** without negative controls for documented butterfly visits and taxa visibility/genus. The signed source is photographic presence, not a surveyed flora or insect absence frame. To establish a biological effect need true effective day/night service, pollen limitation/seed production or independent longitudinal moth monitoring, and robust ground/landcover validation.

## Execution/reproducibility

On the feature branch, the new workflow tests original raster geometry with synthetic pixels, restores the already-verified site-source artifact and authentic high-resolution GSHHG from existing Actions, MD5-checks one official 2024 PCNL raster, and saves actual 2025 photo-site spatial measures as a new non-destructive Actions artifact.

No RGB photo downloads, ROI reclassification, unknown-colour imputation, refitting of old flower-white coefficients, or changes to existing ALAN main results. A successful source download and connectivity table would be a **geospatial exposure result only**, not a confirmed light effect on a population or a nocturnal-pollination mechanism.
