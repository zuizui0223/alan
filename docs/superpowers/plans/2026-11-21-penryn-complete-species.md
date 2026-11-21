# Penryn complete-species analysis implementation plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task by task.

**Goal:** Estimate named-species composition and trophic effect modification from the same 54 physical units.

**Architecture:** Validate complete species against survival and community allocation; build 24 species lnRRs and their empirical joint covariance; apply explicit linear contrasts; render and review results.

**Tech Stack:** Python, numpy, pandas, scipy, matplotlib.

**Spec:** `docs/superpowers/specs/2026-11-21-penryn-complete-species-design.md`

## Global constraints

Use 54 mesocosms, not 216 species rows as independent lighting replicates.
Do not impute missing biomass, change existing outputs or infer a wetness regression.
No raw-data publication changes or main merge.

## Review focus

- Source survival conflicts must exclude a species even if its biomass table is complete.
- Pair species by mesocosm ID before estimating covariance.
- Shared dark controls must create covariance between light-colour contrasts.
- Trophic modifiers require differences of light responses, not comparisons of significance.
- Arithmetic-mean and mean-log estimands must remain explicitly distinct.

## Task 1: validated analysis and results

- [x] Write tests for source selection, duplicate/missing inputs, hand-computed shared covariance and scale invariance.
- [x] Run tests and observe failure because the new analysis module is absent.
- [x] Implement `database/alan_experiments/plant_sources/analyze_complete_species.py` with CLI `--out`.
- [x] Run targeted tests and produce tables, full covariance and source hashes in `plant_sources/complete_species_trophic/`.
- [x] Independently verify original PDF values and effect/covariance arithmetic.

## Task 2: scientific presentation and integration

- [x] Render a multi-panel figure and visually check saved PNG/PDF.
- [x] Write Japanese results and limitations, link from current evidence.
- [x] Obtain independent scientific/code review; resolve concrete material findings.
- [ ] Verify source hashes, outputs and diff; commit and open a draft PR.

## Execution record

Baseline: current main acquired anew; three existing meta-analyses independently
recomputed without using their synthesis functions and matched stored values.
Three existing experiment-DB tests passed. A fresh clone and dedicated analysis
branch provide isolation. The user already authorized continuation; no new
approval gate is needed for analysis and a draft PR.

Task 1 complete: 6 new and 3 existing tests passed; original PDF216 rows, all values and treatment labels match independent Poppler parsing. Source-number strings were preserved through the subset CSV; repeat verification has zero lexical differences.

Task 2 review: independent arm-mean/Jacobian calculations match24 effects and576 covariance entries within1e-14. PDF/PNG inspected after margin corrections; all labels visible. Nominal white upper bound0.996925 is explicitly qualified by the two-colour interval including1. No new pooled effect or independent programme is counted.

Final review: no material findings. Overview wording clarified as the geometric mean of comparator responses. All6 new and3 existing tests pass. Portable independent covariance/source verifiers and their receipts are retained. No unresolved reviewer findings.
