# ALAN A | Source sufficiency and targeted original-data requests (2026-11-15)

## What was actually verified

This is a **post hoc source audit**, undertaken after existing positive/negative published ALAN results were seen. It does not claim a preregistered analysis or source-data meta regression.

- 12 registered plant reproductive candidate publications (9 direct assigned-light experiments, 3 observational light gradients), plus two mechanism-only publications kept outside the reproductive denominator. **Two** of the nine published direct ALAN experiments explicitly crossed artificial light with full/day/night/open pollinator access: Macgregor et al. 2019 and Lu et al. 2026 JPE.
- The confirmation-ready number with a common open-reproduction light/dark log response ratio, true clustered standard error, independently measured same-population dark day/open functional nighttime dependence, and the shared-control covariance is **0**.
- A verified 7-test source-readiness census ran successfully: [GitHub Actions #43703217786](https://github.com/zuizui0223/alan/actions/runs/43703217786), artifact #14156653591, ZIP digest sha256:1d1b5d707fb5c5f3b40a76dd4714a4a2c583342a62350750c31e243a898972c3. Model-ready classifications are in results/alan_A_endpoint_independence_readiness_20261115.json and .csv; logic in scripts/analysis/alan_A_strict_readiness_20261115.py.
- Different public alternatives for the publishers' supplements were actually checked with [GitHub Actions #43702749632](https://github.com/zuizui0223/alan/actions/runs/43702749632), 3 unit tests passed. Lu accepted manuscript / Silverchair static DOCX / Macgregor Wiley exact Appendix S1 remained HTTP403; Chinese journal metadata HTML200, without underlying arm tables. Zero new original raw treatment records were acquired. The workflow records this access fact, NOT a negative biological finding.

## Original studies closest to factorial inference

| Original paper | Known true design | Published result with proper treatment scope | Missing data that block beta_V |
|---|---|---|---|
| Macgregor et al. (2019), Ecosphere, DOI 10.1002/ecs2.2550 | 5 field locations × 4 temporal rotations, 5 light levels × four access regimes; 169 female plants, 417 flowers | Full-night illumination increased successful pollination; **original light×caging interaction was not statistically significant**, not proven precisely zero. Published day=37.6%, night=52.8%, open=72% are **pooled across all lighting levels**, not dark-only D0/O0 | Source scripts refer to Data/SeedSetYN.csv, Data/SeedSetBinom.csv, Data/SeedCount.csv, none included in the six-file GitHub source tree or inspected Zenodo ZIP. Need site/round/access/lamp arm counts and covariance. |
| Lu et al. (2026), Journal of Plant Ecology, DOI 10.1093/jpe/rtag160 | Five dark/lit plot pairs, two plant species, 2×4 light/access factorial | **−56% P. kansuensis developed seeds/ovules is NIGHT-ONLY light vs dark**, not open total reproductive yield. Open seed effect is nonsignificant with unknown point estimate. For same pollen outcome, open lit/dark 0.7691, night-only lit/dark 0.4620; ratio-of-response-ratios 1.665, log +0.510, without SE | All 16 light×access×plant original reproductive cells, source Supplementary Table S4b and plot-level clustered covariance. Outcome is developed seeds per ovule in a formed fruit, not per original assigned flower. |

**Even if original data for BOTH factorial studies are recovered, they comprise only TWO independent paper-level clusters**, insufficient by themselves for a precise and generalizable global V slope. The remaining seven direct ALAN studies would also require same-population, independent functional day/open dependence rather than assigning a default V by flower morphology or genus. The number of published effects is not the number of independent illuminated sites.

## Minimum data request specification: NOT SENT

**For Lu et al. 2026**: 5 light/dark field pairs × 2 plants × four access regimes (full exclusion, day only, night only, open) with the true lighting-site pairing, plant and date IDs, originally assigned flowers, mature fruits, developed seeds, total ovules, model link and random effect, light×access coefficients and standard errors. If plant-level source rows cannot be shared, a table of per independent plot-pair × access × plant success/trial counts, with source covariance and published model output, could provide reproducible arm contrasts. The original accepted manuscript's Supplementary Table S4b is especially important.

**For Macgregor et al. 2019**: the three named original input CSVs in the authors' six-R-script repository, or independently aggregatable treatment cells with the 5 lighting types × four caging treatments × four experimental rotations × site ID. Include seed-positive successes/total flowers, plant reuse identifiers, share of consumed capsules vs usable seed counts, factorial GLMM output, site and run covariance, and Appendix S1 Tables S3–S7.

The requests are a research work product; **no email to authors was sent**.

## Causal interpretation limits that cannot be relaxed merely to gain k

- Lu's 35.71% nightly contribution and 56% night-only illumination seed decrease are **not on the same unrestricted reproductive contrast**, so the previously proposed extra 20.29 percentage-point daytime loss was invalid and remains retracted.
- Grognuz & Knop (2026) demonstrates an ALAN-associated seed change after preceding nights when target flowers were bagged, then given only a daytime pollinator opportunity: an alternative daylight carryover pathway, not a study with confirmed nocturnal guild depletion causing all seeds.
- Macgregor and Li report positive reproductive directions under some light configurations; do not select only negative examples. Boom's light-sensitive nursery moth changed fruit infestation and ovule fertilization but the author's estimated final seeds/flower ANOVA was p=.63; do not call a significant net loss.
- The main hypothesis survives as **effect modification of total plant reproduction by genuinely independently measured nocturnal functional dependence**, not the exclusive necessary cause of all negative ALAN impacts. The island literature describes channel composition; there is no demonstrated island×ALAN reproductive experiment interaction in this source set.

## Disposition

**Continue the strict primary A meta-regression only after reproducible original arm means, the lighting assignment cluster covariance, an independent V for the same population, and adequate independent replication are available.** Until then, preserve this as an evidence map and a source-request programme. The coefficient beta_V, its 95% CI, and the predictive island effects remain unestimated.
