# Murphy2026 grass competition design review

Reviewed 2026-11-20. Independently read full Methods and Results from the [author-hosted published PDF](https://smurphlab.weebly.com/uploads/1/2/6/5/126523979/murphy_et_al._2026.pdf), [DOI 10.1016/j.biocon.2025.111586](https://doi.org/10.1016/j.biocon.2025.111586). Raw-data location supplied for this study: [Zenodo 17479452](https://zenodo.org/records/17479452); raw observations were subsequently inspected as documented below.

## Primary-source summary (under 200 words)

The 2015–2016 greenhouse experiment crossed two native grasses, two exotic grasses, three exotic densities and two lighting conditions, using 240 conetainers. One greenhouse room received ALAN from six sodium lamps; the other received ambient light. Authors explicitly acknowledge absent independent replication of lighting. Flats and pots rotated within treatments, not between rooms. Harvest followed thirteen weeks.

At 1:1 planting density, exotic biomass responded more strongly to ALAN than native biomass. However, ALAN did not significantly modify the negative relationship between native composite growth and exotic composite growth. It did modify the native C:N relationship with exotic growth. These are different endpoints. The composite growth measure combined culm number, length and dry mass; native and exotic PCAs were fitted separately. No belowground biomass or herbivore response was measured.

## Independent eligibility assessment

The greenhouse setting itself is allowed by the current inclusion protocol. The limiting issue is **one independently illuminated room per arm**, not an indoor-behaviour exclusion. Six lamps are components of a single room-level treatment; ten pots per species-density-light combination, flats, individual culms and repeated survival observations are subsamples. Random placement of competition combinations and rotation within a room do not replicate lighting or separate its effect from room conditions. Similar recorded temperature/humidity do not eliminate unmeasured room confounding.

Consequently, a conventional treatment-arm sampling variance based on pots would quantify subsampling uncertainty while omitting variation among independent light environments. It cannot supply the required independently replicated ALAN effect for ordinary inverse-variance pooling. Keep this study in a clearly labelled unreplicated experimental/contextual stratum; raw data recovery cannot fix this design limitation. Do not set n=6 from lamp count, n=10 from factorial-cell pots, or infer an independent paired-light design.

## Distinguish the reported estimands

1. The 1:1 density comparison concerns aboveground biomass differences by nativity and lighting. Native and exotic plants sharing a pot are dependent. This is not evidence that native biomass declined under ALAN, nor a comparison of competitors present versus absent.
2. The main competition models use realized exotic growth PC1, a post-treatment response, rather than solely the assigned density. A light-by-realized-growth association does not independently establish a causal mechanism. Conditioning on realized growth and centering native responses within lighting treatments further separates this estimand from raw biomass lnRR.
3. There is no zero-competitor treatment: assigned densities are one, two or four exotics per native. A presence/absence competition contrast or its response-ratio interaction therefore cannot be reconstructed without inventing a reference.
4. Native C:N, biomass and leaf toughness are distinct responses. Changed C:N can indicate altered tissue chemistry; it cannot be counted as measured herbivore consumption, plant mortality or reproductive loss. The native-growth null modification must accompany the C:N finding rather than disappearing behind the abstract's general competition claim.
5. A coefficient for a PCA score or light-centered response is not an arm mean. Native and exotic PCA axes have separately estimated scales. Do not exponentiate these coefficients as biomass ratios or infer arm SDs from published F statistics.

## Decision

Design eligibility for ecological setting: yes. Eligibility for the replicated-light quantitative synthesis: no, owing to unreplicated room-level allocation. Preserve the source and its competition results with this restriction, but add no ordinary independent-unit effect rows. Author-model findings may be described separately, without calling them an independent test of the macroanalysis wet-affinity mechanism. Check shared greenhouse/programme provenance with Murphy2022 before combining contextual reports; this review did not establish raw-data overlap.

## Independent original-workbook identifier check

Inspected `Murphy et al. data 2025 ALAN competition.xlsx` directly with openpyxl/pandas. Retained Plant-data rows with nonmissing `Species_TRT` rather than counting formatted worksheet extent. Normalized `Rep` to lowercase for cross-sheet comparison and used `(Species_TRT, Density, Rep, Trt)` as the pot key, adding `Status` for chemistry.

- Plant data: 240 rows and 240 unique pots; each of 24 species-pair × density × light cells has ten pots.
- C_N data: 397 rows but 379 distinct pot/status combinations. Every chemistry pot maps to the plant table.
- There are 36 chemistry rows in 18 duplicate-key groups. Four groups have identical numeric `%N`, `%C`, and `C:N`; fourteen groups have two different values for all three variables.
- Relative to the nominal 240 pots × two nativity statuses, 101 pot/status combinations are unrepresented. This counts absent rows, not necessarily missing laboratory observations from extant plants; causes are unresolved.

For example, BGCG/density1/rep h/NL/NATIVE has C:N values 21.408067 and 15.094896. BGSB/density1/rep i/L/NATIVE repeats the same C:N value 29.987426. These contrasting duplicate types rule out assuming all duplicates are either accidental copies or technical replicates. Preserve original rows and resolve the sample-ID/replicate provenance before merging chemistry into a one-row-per-pot/status analysis. No rows were deleted, averaged or reassigned here. The original source files were unchanged. These checks corroborate the parent's reported counts, while the unreplicated-room limitation remains decisive for ordinary ALAN meta-analysis.
