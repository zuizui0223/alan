# Full-text batch 4 — 2026-11-20

Fifteen newly retrieved primary full texts were individually reviewed using Methods, Results and original-data availability/supplement statements. The batch comprises the requested fourteen IDs plus late-retrieved PMC3861544. Same column schema as batch 1.

- 3 eligible designs, extraction/verification status tracked separately.
- 9 held for independence or scope.
- 2 observational-light contrasts, separated.
- 1 indoor behavior-only study, separated.

**No new effect estimates were extracted in this screening batch.** Existing Parkinson2020 and Giavi2021 extractions remain eligible search records and are not discarded as duplicate articles merely because extraction already exists.

## Useful additions and boundaries

**Field firefly mating, PMC9364009, doi:10.1098/rsos.220468:** actively dialoguing pairs received a 5-lux flashlight treatment or dark sham; 34 field pairs over 14 evenings. Original CSV/code are in the supplementary ZIP. Prioritize the actual field mating endpoint, verifying arm counts, unique animals and session clustering. The article also includes binary-choice and laboratory assays: these are different experiments, not extra independent copies of the field result. The outcome is binary, so continuous-response lnRR extraction is not automatically appropriate.

**Plant performance, PMC10613542, doi:10.1098/rstb.2022.0358:** [Dryad hhmgqnknt](https://doi.org/10.5061/dryad.hhmgqnknt) contains original data. Twelve indoor EcoUnits represent ten ALAN settings with two repeated extremes and single intermediate doses. The first harvest declines with light; the second harvest shows little change. This is genuine ecological biomass and diversity evidence, not behavior-only laboratory work. It remains an indoor-scope/replication adjudication case, and shares the Cesarz Ecotron 2020 programme and EcoUnits; do not count it as a new independent programme.

**Plant disease, PMC3861544, doi:10.1371/journal.pone.0081282:** night illumination suppresses basil downy-mildew sporulation and net-house epidemics. This adds an important antagonistic interaction: light can reduce pathogen damage rather than increase damage to plants. Figure 10 supplies means/SD and n=13 containers, but those containers occupy one contiguous lit half-row and one dark half-row in a single net-house; sequential trials share the inoculum system. The Methods and Figure 10 reverse which compass-side half was illuminated. Hold independent-treatment extraction pending allocation clarification. Growth-chamber detached-leaf assays must be separate from net-house ecological outcomes. This was not excluded as indoor behavior.

## Existing eligible studies confirmed

- **PMC7544032, Parkinson2020:** twenty independent littoral mesocosms, ten with their own added lights and ten without. Original [GitHub data](https://github.com/LizzParkinson/ALAN_Data), already extracted in the repository. Fish manipulation failed and was dropped; do not interpret this as a successful fish-predation factorial. Spider abundance, mass and web-prey counts are distinct endpoints, and web captures do not supply an attack rate. Phone-based reported light measurements require caution.
- **PMC7966740, Giavi2021:** six paired lit/dark meadow blocks with repeated matched-effort daytime surveys; [Zenodo 4540407](https://doi.org/10.5281/zenodo.4540407), already extracted. Responses differ among selected plant species and visitor orders. Counts do not establish pollen-delivery compensation. Check programme/year identities before combining with other Swiss reports.

## Independence and exposure checks that changed interpretation

- **PMC11143364:** Westhavelland has one illuminated field and one dark field, each with twelve lamp/trap subsamples; municipal retrofit arms are all lit. Original [Figshare 23773365](https://doi.org/10.6084/m9.figshare.23773365). Shares the Westhavelland programme, and municipal data overlap the parasitoid retrofit report in batch 1. Monthly catches and changing lights across years cannot create independent treatment fields.
- **PMC10997633:** one boat light in one lake, repeatedly switched on/off over five nights. Vertical redistribution of *Chaoborus* is evident, but total abundance did not differ by lighting condition. Independent installation variance is absent. Separate anoxia-survival experiments do not measure mortality caused by ALAN.
- **PMC3961812:** eight large indoor ecological stream flumes on a lighting gradient, one per dose and one covered dark control. Not randomized across positions; weeks substitute for spatial replication in the article. Retain null drift findings and the indoor ecological scope flag; do not classify this as behavior-only.
- **PMC10265009:** field marine-isopod arenas compare direct/diffuse light and shadows; all quantitative arms are illuminated. Camouflage and choice behavior do not measure predation rates.
- **PMC4376897:** bridge-mayfly field observations and polarized/unpolarized torch experiments. Quantitative arms are lit; brief off periods reset the swarm rather than provide a sampled dark treatment. Important mechanistic reflection evidence, held separately from replicated ALAN/dark ecology.

## Other scope decisions

**PMC8187998** transplants fish/anemones into existing illuminated versus dark reef sites; **PMC8397135** releases or monitors females in existing road-light environments. Manipulating organisms does not randomize or manipulate the lighting contrast. Both are preserved in the observational-light layer, with original repositories recorded. These decisions are about exposure design, not the taxa.

**PMC12041495** tests oyster valve rhythms in two tanks through sequential light schedules: indoor behavior-only. **PMC12339407** includes fly size/survival and molecular traits under laboratory lighting: held for indoor fitness scope and unverified independent light units, not automatically behavior-only. **PMC12734188** combines natural-skyglow field orientation with indoor gradient manipulations; it does not provide a manipulated outdoor dark-control population response.

No master tables changed and no commit. CSV records contain the source sections, replication qualifications and overlap flags for all fifteen papers.
