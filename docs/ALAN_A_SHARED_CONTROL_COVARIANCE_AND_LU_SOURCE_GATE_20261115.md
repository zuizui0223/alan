# ALAN A: shared-control bias and Lu 2026 raw treatment source audit (2026-11-15)

## Status

This audit was developed after seeing positive and negative original ALAN outcomes. It is a post hoc method and data-availability audit, NOT retrospective preregistration. The ALAN night-dependence meta-regression coefficient remains **not estimated** because its joint moderator, reproduction endpoint and variance are not yet source-complete.

## 1. A negative slope can arise from reusing the same dark/open denominator

Define V = 1 - D/O, where D is dark with day-only visitors and O is dark with day+night visitors. Define y = ln(L/O) where L is lit with day+night visitors.

By first-order delta approximation assuming independent treatment arm sampling:

Cov(V,y) ≈ (D/O²)(-1/O) Var(O) = -D Var(O)/O³ < 0 when D>0.

Even when ALAN truly has zero reproductive effect on every plant, sampling noise of O pushes V and y in opposite directions. A naive meta-regression can thus manufacture the preselected negative beta.

### Real numerical simulation execution (synthetic data, NOT ALAN evidence)

Source: GitHub Actions run #43694065051 SUCCESS: https://github.com/zuizui0223/alan/actions/runs/43694065051 . Test count 4, all passed, pinned numpy 2.1.3. Artifact #14153651325, ZIP sha256 0dbab68bdd94f42fa336190433772e365ae617b473c6f50f3df283ec3520f3ad. Reproducible source: scripts/analysis/alan_A_shared_control_null_20261115.py .

Synthetic assumption: eight virtual plant species, true V between 0.15 and 0.75; dark and lit open probabilities both 0.60; zero actual light effect; 50,000 simulated datasets for n=20,60,200 each.

| Flowers per treatment arm | Shared-O mean estimated beta | Fraction beta negative | Independent-O mean beta | Independent-O fraction negative |
|---|---:|---:|---:|---:|
| 20 | -0.1988 | 0.7305 | -0.0017 | 0.5012 |
| 60 | -0.0844 | 0.6556 | +0.0015 | 0.4974 |
| 200 | -0.0286 | 0.5939 | -0.0003 | 0.5006 |

This is a **statistical diagnostic** of the proposed estimand, not a result about actual ALAN or proof of any particular biological negative coefficient. Size depends on true variance in V, baseline reproductive success and actual n. Split samples are independent only if split at independent experimental units, not arbitrarily split flowers on one light-treated site.

### Required inferential correction

1. When a paper reports V and y using the same dark-open observations, record the covariance, and fit measurement-error-aware multivariate hierarchical synthesis rather than naive regressions.
2. If V is independently measured on a separate site/date, record the biological external-validity mismatch and moderator uncertainty; do not pretend a different population has exactly the same V.
3. True units are light-assignment clusters, blocks, reused plants and periods; flower counts are not independent lamps.

## 2. Lu et al. (2026) J Plant Ecol: public source-access execution

DOI: 10.1093/jpe/rtag160. Source-verification run #43694005508 SUCCESS (three unit tests): https://github.com/zuizui0223/alan/actions/runs/43694005508 . This success certifies checking public endpoints, NOT accessing the raw experiment. Probe source: scripts/analysis/alan_A_lu2026_source_probe_20261115.py ; verified artifact #14154070962, ZIP sha256 199f4d66138f104d38939cf25664c8c7caade356c69d6e7188fe72846c2ae896.

- JPE official journal public article page: HTTP 200, 70,022 bytes, but no directly downloadable individual treatment-data link.
- Oxford Academic publisher landing and two standard supplementary DOCX paths: HTTP 403.
- Oxford Academic accepted PDF path: HTTP 403. Crossref metadata HTTP 200 exposes official PDF location with 68754751 in URL but no accessible quantitative data from this check.
- OpenAlex metadata HTTP 200, is_oa=true but best_oa_location has no pdf_url.
- Europe PMC metadata HTTP 200, no PMCID in this record.

Therefore **all 16 reproductive treatment cells remain unacquired** (2 plants × 2 lighting × 4 access conditions); numerical interaction beta/SE, source-unit n and shared-control covariance are NOT recovered. **Important later source-scope correction:** IWE's independently audited accepted manuscript Results p14/Discussion p16 shows the author-reported -56.0% Pedicularis developed seeds/ovule difference belongs to **night-only access**, NOT open day+night reproduction. The open seed change was reported as not significant with numerical effect UNKNOWN, not zero. The -56% and 35.71% night-contribution fraction refer to different access contrasts and cannot be subtracted to infer additional daytime loss. These remain author results, not a matched two-species inferential interaction. Verified 12-row IWE fixed-source crosswalk: https://github.com/zuizui0223/alan/blob/analysis/alan-floral-timing-dependence-20261115/docs/ALAN_A_LU2026_IWE_SCOPED_NIGHTONLY_CORRECTION_20261115.md . An HTTP 403 is not absence of research and is not a biological negative result. No access-control bypass attempted.

### Original-author data request specification (not sent)

Request 16 cells of seed/ovule and fruit/flower counts, n of independently assigned lighting installations and matched blocks, original plot/plant/time IDs, sample/SE or residual covariance, and exact seed-reproduction lighting×day-night-access coefficient with link scale and SE. Aggregates by genuinely independently assigned lighting unit are a useful privacy-preserving alternative to raw coordinates.

## 3. Diagnostic direct experiment contrast, distinct from the main cross-study slope

When source factorial cells exist, define within-study difference-in-differences:

Delta_within = ln(R_lit_open/R_dark_open) - ln(R_lit_day/R_dark_day).

Negative Delta is consistent with a larger ALAN reproductive loss when nocturnal access remains available than when day-only access is imposed. It is NOT a proven moth-mediated causal mechanism if direct plant light responses, fruit herbivory, timing or daytime visitors change too. This is a separate mechanistic diagnostic, not a replacement for cross-study beta_V.

The original Macgregor et al. 2019 study reports the lighting×caging pollination-success interaction as non-significant, but its numerical interaction coefficient and authentic cell-level standard error remain unrecovered. Lu 2026 has the factorial design, but this audit cannot extract its full reproductive cell table.

**Next gate:** genuine ALAN lighting-by-access treatment data with source n and uncertainty; do not infer an empirical beta_V, island×ALAN interaction or nighttime-moth light susceptibility from abstract percentages. This entire source audit is written on Draft PR #78 only, main unchanged.
