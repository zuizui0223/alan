# ALAN A | Lu 2026 × IWE cross-repo verified treatment scope (2026-11-15)

## Major corrected scientific conclusion

Lu et al. (2026), Journal of Plant Ecology, DOI [10.1093/jpe/rtag160](https://doi.org/10.1093/jpe/rtag160), manipulated nocturnal lighting in **five pairs of lit/dark plots** crossed with **four pollinator-access treatments** in two alpine plant species. The independently audited 33-page accepted-manuscript Results (p14) and Discussion (p16), stored in the user's IWE repository, establish:

**The published -56.0% Pedicularis kansuensis seed-set decrease pertains to the NIGHT-ONLY pollinator-access arm.** It is NOT the ALAN effect on natural, unrestricted day-and-night plant reproduction. Seed set in this study was the proportion of developed seeds among ovules in a formed fruit, not mature seeds per original flower.

The open-access seed-set lighting response was reported by the authors as **not statistically significant**. Its numerical effect estimate and standard error were not extracted; "not significant" does NOT imply effect exactly zero, complete compensation, equivalence, or an effect with adequate power to rule out harm.

### Quantitative original text values, correctly matched

| P. kansuensis, same endpoint | ALAN vs dark relative change | Light/dark ratio | Natural log of light/dark ratio |
|---|---:|---:|---:|
| Conspecific pollen per stigma, open (day+night) | **−23.09%** | **0.7691** | −0.26258 |
| Conspecific pollen per stigma, night-only | **−53.80%** | **0.4620** | −0.77219 |
| Developed seeds per ovule in formed fruit, night-only | **−56.0%** | **0.4400** | −0.82098 |
| Developed seeds per ovule in formed fruit, open | **not significant** | **unknown** | **unknown** |

The point-value access-specific pollen response-contrast is:

B_pollen = ln[(Y_lit,open/Y_dark,open)/(Y_lit,night/Y_dark,night)]
         = ln(0.7691/0.4620)
         = **+0.509656** (response-ratio ratio **1.664719**).

This is arithmetic derived from **rounded published percentages**, NOT a new GLMM interaction coefficient, independent source pot-level fit or causal mediation percentage. In particular 1.665 does not mean extra absolute pollen is 66.5%, nor that daytime pollinators supplied 66.5% of pollen or maintained seed output. It means the *reported relative light response* of open pollination was less negative than night-only pollen receipt.

The 2 plants × 2 endpoints × 3 access conditions yield **12 original source-claim slots**. Three numerical light/dark ratios are explicitly reported, nine comparisons described as non-significant have *unknown* point effects. For the four possible same-species×same-endpoint open-versus-night paired contrasts, only the P. kansuensis pollen comparison has a numeric point estimate. **No seed open/night effect modification is numerically estimable; no source sampling covariance is available**. The original 16-cell seed experiment (two lights × four access arms × two species) and associated plot-level outcome variance remain inaccessible in the present source probe.

## Retraction of the erroneous 35.71% versus 56% arithmetic

Previously, the 35.71% relative night-pollination contribution to seed set was compared with "56% ALAN loss" under an assumed common unrestricted reproductive denominator, yielding a hypothetical extra 20.29-point daytime loss. **That inference is invalid** because the 56% pertains to the **night-only access arm** and is a conditional developed-seed/ovule fraction, while 35.71% refers to a published nocturnal share using another source contrast. No 20.29-point excess loss or 36.23% mediation estimate survives this access-scope correction. The code returns null for that hypothetical bound in the original P1/P3 report.

## Why this is directly relevant to the main A hypothesis

- The original main A target is **lit/open versus dark/open reproduction**, as a function of independently measured marginal night dependence V. Lu's -56% night-only observation **cannot serve as the main y = lnRR(lit/open, dark/open)**.
- Lu is evidence that ALAN impairs the isolated nocturnal channel and a suggestive point-level example that unrestricted pollen reception can buffer that impairment. **It is not yet a statistically supported ALAN×night-dependence interaction on final seed yield**.
- Macgregor et al. 2019 is another factorial lighting × day/night access study but did not find a statistically significant interaction, and its published 72% open, 37.6% day, 52.8% night rates are pooled across lighting (they cannot estimate dark-only V).
- The two experiments together provide a real *experimental mechanistic roadmap*, not a globally generalizable slope; their original outcome denominators (developed seed/ovule vs any seed/flower) differ and source covariance is not recovered.

## Independently reproducible cross-repository verification

Source IWE user repository: [lu2026_access_scoped_contrasts_v01.csv](https://github.com/zuizui0223/iwe/blob/main/data/registry/lu2026_access_scoped_contrasts_v01.csv), fixed commit `cd0130b4730735f0ec90b1162780ad67b19a222a`, Git blob `73a04ebd962a3113ebfae410c1720468094ab60a`; original transcription file SHA256 `b8f18e68f3efcbdc108b59fcd79df9551fd6651cea6a34a0f6ae9fdf5d7b6b43`. Source textual page map: [IWE fulltext scope](https://github.com/zuizui0223/iwe/blob/main/data/registry/lu2026_fulltext_scope_source_v01.json). The 33-page publisher PDF could be parsed from a cached web copy in IWE, but the PDF binary, its original numerical Table S4b and PDF screenshot were not retrieved; no binary PDF hash is claimed.

New ALAN analysis script: `scripts/analysis/alan_Lu2026_iwe_access_scope_sync_20261115.py` (fetches exact pinned IWE audit and checks source SHA, 12 row keys, access, endpoint, absent SE and 3 published numeric ratios); unit tests: `tests/test_alan_Lu2026_iwe_access_scope_sync_20261115.py`; workflow: `.github/workflows/alan-A-Lu2026-IWE-access-scope-source-20261115.yml`.

[GitHub Actions run #43701694874 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43701694874): **7 independent source-guard unit tests passed**, source SHA check succeeded and gave 12 slots / 3 numerical light ratios / 1 numerical pollen open-v-night point / **0 numerical seed open-v-night pairs** / 0 A beta_V. Source-only report artifact #14156457362, SHA256 `17ec60cd0f73169fadfc04593ee97d6b1384fd9cb26cb9c750fdc4fc17273280`.

Scientific summary: *the direction of isolated nocturnal pollination harm is well documented in Lu 2026, but a decline in unrestricted whole-plant seed formation is not quantitatively established by this -56% statistic*. The A meta-regression β_V and island × light effect remain unestimated.
