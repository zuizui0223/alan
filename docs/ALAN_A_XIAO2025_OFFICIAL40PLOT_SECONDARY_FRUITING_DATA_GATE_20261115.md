# ALAN A secondary reproductive phenology and true original supplement audit: Xiao et al. 2025 (2026-11-15)

## Verified publication and source

**Xiao, Wang, Oduor, Wang, Zhang & Liu (2025)**, Biological Conservation 311, 111472; DOI [10.1016/j.biocon.2025.111472](https://doi.org/10.1016/j.biocon.2025.111472). Native Inner Mongolia grassland original manipulation using **40 field plots**, each 4 m², in an **ALAN (yes/no) × animal-exclusion (yes/no)** design. Experimental light was directly manipulated.

**Do NOT relabel animal exclusion as day-only vs night-only effective pollinator exclusion.** The excerpted publisher Methods do not identify the excluded animal taxa and whether pollinators entered. The true independent lighting-installation/block assignment number and plot distribution across lighting installations remain unconfirmed. Forty plots are not automatically forty independent lamps.

The study measured an annual herb, *Elsholtzia densa*. The original Results/Abstract state ALAN advanced the onset of budding, flowering, fruiting and seed maturation by **3.4–6.2 days**, reduced the **number of fruiting inflorescences**, increased specific leaf area by **9.2%**, and decreased secondary branching by **44.3%**. These latter changes are original published results, not effects derived in the new reproduction meta-analysis. The plant trait and fruiting phenotype show a possible direct physiological/temporal ALAN pathway; mediation through moth pollination, predator exposure or herbivory was not isolated here.

**Outcome distinction:** a plant-level number of inflorescences bearing fruits is NOT the probability that an initially marked flower forms a mature, seeded fruit. Nor is an advanced fruiting date an actual viable seed yield change. This original paper is a legitimate direct manipulated ALAN study with a **secondary plant reproduction/phenology endpoint**, but it lacks the fixed A primary flower-level binary fruit outcome and the independent nighttime insect dependence moderator V. Preserve it in the evidence map and in stage-specific future analyses, not in the primary y regression.

## Genuine source retrieval done, not just a search

[GitHub Actions #43708104300](https://github.com/zuizui0223/alan/actions/runs/43708104300) successfully checked the official Elsevier supplementary endpoints for this DOI and the Oxford publisher/Cal Poly Wilson 2021 materials. The official supplemental **DOCX** is actually downloadable:

- URL: https://ars.els-cdn.com/content/image/1-s2.0-S0006320725005099-mmc1.docx
- HTTP **200**, exact file size **315,350 bytes**.
- Original file SHA256 **127950cecd9ac78c0927b89df06b307e92a7a9445b0c054aa40a2e9d272e8a9f**.
- Parsed the real zipped Word document structure, not its filename: **5 paragraphs, 0 Word tables and 3 embedded figures**. The few figure captions describe, among other things, the spectral composition of artificial lighting and top-inflorescence morphology. **There are no original seed/fruit treatment-cell values or model effect/SE tables in this supplement**.
- No source figures/images or original DOCX are redistributed into the repository. No numerical response extracted from figures was invented.
- Tests: [Xiao official DOCX table source CI #43708383911](https://github.com/zuizui0223/alan/actions/runs/43708383911), original SHA pinned, 3 structural parser tests passed. Artifact #14159671136, ZIP SHA256 d58cfa360d6194ac104e5b45e0f2da32bf87b83606b2a772b4c222289827c8f7.
- Follow-up CI also checks that the source registry treats this phenotype as an explicitly secondary, *not* an already extractable first-flower fruit set outcome, and does not grant 40 independent lamp replicates or a nocturnal insect access V.

**Wilson et al. 2021**, DOI 10.1093/icb/icab010, remains a genuinely randomized experimental yucca fruit source (five sites, 15 lit/dark transects, 58 analyzed plants), but the official Oxford supplement ZIP remains HTTP403. Cal Poly thesis landing is accessible HTTP200; no additional original randomized transect fruit data were recovered from it in this probe. This does not make the existing source response zero.

## Main A implications and source-request gate

The published Xiao study supports a competing pathway: ALAN can change plant flowering and reproductive resource allocation, with potential consequences for fruiting independent of a directly measured moth-loss process. It does **not** prove a causal plant-physiology-mediated reduction in final plant fitness, because direct plant and interaction pathways were not quantitatively partitioned in accessible data.

For numeric analysis, request the original **40 plot treatment table** with plot/light tower/block IDs, actual identity of excluded animal groups, per-plant starting inflorescence and flower counts, final fruiting inflorescences, developed seeds and ovules, dates, exact light exposure, model link, coefficient and clustered SE. Without these, no ALAN fruit/seed log response ratio with valid variance or test of ALAN×night-insect dependence can be reconstructed.

This is an evidence **access and design** milestone, not a new independent source study count. The current A register remains 12 reproductive candidates (9 with directly manipulated ALAN, 3 observational gradients), none of which currently provides jointly source-verified open-plant reproductive logRR and same-population dark-only V with covariance. No global A beta_V has been estimated.
