# Knop et al. 2017: original-source extraction audit

Audit date: 2026-11-20. DOI: [10.1038/nature23288](https://www.nature.com/articles/nature23288).

## Status

Two original, meadow-level visitation summaries extracted; paired covariance and raw site records remain unavailable. No pool-ready paired variance and no independent-meadow fruit effect have been reconstructed. This is not a claim that the raw data do not exist. The article explicitly offers data and code on request; no public raw deposit was verified in this search. No author contact was made.

## Sources and design

The original article, Methods and Extended Data Tables were read from the unabridged Nature paper reproduced on pages 172–187 (1-based) of a [government-hosted compilation](https://planesynormas.mma.gob.cl/archivos/2021/proyectos/ANTECEDENTES_ENTREGADOS_POR_ROC-OIKONOS.pdf). Original [Supplementary Tables 1–3](https://static-content.springer.com/esm/art%3A10.1038%2Fnature23288/MediaObjects/41586_2017_BFnature23288_MOESM5_ESM.pdf) were separately obtained from the publisher. Source checksums are in `source_manifest.json`; PDFs and rendered Extended Data Tables are in the external knop2017 cache.

The 2014 study is an eight-meadow unlit day–night network, not an illumination experiment. The ALAN experiment is June–September 2015: 14 independent ruderal meadows, seven illuminated and seven dark, arranged as seven proximity pairs and sampled simultaneously within pairs. Light was 4000 K LED on 4 m poles, switched on below ambient 5 lux. Mean measured light was 52.0 lux within 5 m, 25.7 lux at 5–10 m and less than 3.9 lux beyond 10 m. Light allocation randomization is not explicitly established by the Methods; random plant bagging is established. Do not conflate these.

For visitation, a 100 m transect was walked every 30 minutes for four hours, eight rounds/night, on five sampling nights, except one illuminated and one dark meadow sampled four times. The published summary pools the observation window. Meadow, not insect, transect pass or sampling night, is the independent lighting unit. The final-time rule cannot be fulfilled from these pooled summaries; preserve the actual time window and do not label this final-time extraction. The supplementary tables contain coordinates, plant names and pooled visitor taxon percentages, not site-level visit counts or reproductive numerators/denominators.

## Actual extracted summaries

| Endpoint | Dark mean ± SEM | Lit mean ± SEM | Independent n per arm | ln(lit/dark) |
|---|---:|---:|---:|---:|
| All nocturnal visits, full ALAN dataset | 105.9 ± 22.3 | 40.1 ± 8.8 | 7 meadows | −0.971119 |
| Cirsium oleraceum nocturnal visits, full ALAN dataset | 74.6 ± 23.3 | 27.2 ± 7.5 | 7 meadows | −1.008924 |

Sources: Extended Data Table 1 and 4, visually inspected from original pages. SD was converted as SEM × sqrt(7), giving dark/lit 59.000254/23.282612 for all visits and 61.646006/19.843135 for Cirsium. `visitation_summary_primary.csv` retains means, SEM, derived SD and provenance. These rounded tabulated values are not raw counts. Neither row is currently approved for the paired meta-analysis: the paired covariance is missing. The stored independence-reference variances 0.092501 and 0.173582 are explicitly NOT verified paired variances. For seven pairs, the correct delta variance additionally subtracts 2×covariance/(7×mean_lit×mean_dark).

All-visits and Cirsium-visits outcomes overlap biologically: Cirsium is a subset of the whole community. Do not count the two rows as independent evidence. Restricted pollen-informed visitor subsets are also overlapping outcomes, not additional experiments. The full visitor metric is flower contact, not necessarily pollen deposition or successful pollination.

## Fruit experiment: denominator and model traps

Five lit and five dark meadows (five meadow pairs) were used. Each meadow had five matched pairs of plants; one plant per pair was randomly bagged before flowering. Two heads per plant were counted. Thus each light × bagging arm has 25 plants and 50 heads, but only five independent meadows. Fig. 4 reports head-level means ± SEM (n = 50), which cannot be relabeled as meadow-level SD or independent n = 50.

Developed fruits and aborted ovules were counted. Fruit success proportion requires the developed/(developed+aborted) denominator and hierarchical binomial or appropriately clustered log-odds treatment. No blind lnRR of a percentage or conversion of the reported 13% loss to a complete effect is justified. Original models nested plant pair within meadow within meadow pair, and included an observation-level random effect for overdispersion.

Extended Data Table 3 reports treatment × bagging interactions. A coefficient for that interaction is NOT the direct lit-versus-dark contrast among open flowers. A valid open-flower contrast requires adding the light main coefficient and interaction, with their covariance to derive its uncertainty. This covariance and meadow-level reproductive totals were not available. The paper's interpretive claim of incomplete compensation by daytime pollinators therefore remains qualitative experimental evidence in the database, not a reconstructed independent numerical effect. The study did not perform separate daytime-only versus nighttime-only exclusion here.

## Relationship to Giavi 2021

Giavi's experiment occurred in 2016, not 2015, and used 12 meadows in six pairs. `sites_2015_primary.csv` transcribes all 14 Knop experimental coordinates from Supplementary Table 1. `giavi_coordinate_crosswalk.csv` matches every Giavi coordinate from the existing primary-coordinate table to its nearest Knop coordinate using haversine distance.

Eight of the 12 Giavi sites are within 0.60–33.93 m of Knop sites, with matching light treatment; four are 15.65–20.38 km from their nearest Knop coordinate. The 50 m reporting threshold is a proximity flag, not proof of identical meadow boundaries. Pair identities differ for most matched sites; only Giavi pair 1 exactly preserves the corresponding Knop pair in this nearest-coordinate matching. This supports probable shared meadow programme across different years rather than duplicate observations. Conservatively cluster Knop and Giavi at programme level and record site/year dependence; do not exclude one as a duplicate, and do not assume fully independent studies. Lighting history for matched Giavi meadows needs explicit verification before coding initial exposure.

## Ecological applicability and remaining needs

Cirsium oleraceum and the listed companion plants can provide wet-affinity biological context, but the experimental Methods describe ruderal meadows, not measured water-table conditions or aquatic margins. Wetland status, water distance and soil moisture must remain unclassified until coordinate-based covariates or primary site measurements resolve them. This experiment supports a pathway from lighting to fewer nocturnal flower contacts and altered reproduction; it does not directly test water reflection, spiders, predator-mediated losses, long-term plant population decline or island effects.

To finish a numerical paired extraction: obtain 2015 visit counts and effort by meadow/date, reproductive developed/aborted counts by meadow/plant/head, and the meadow pairing key. Those data permit effort-standardized visits, paired covariance, final-window analysis if warranted, and a clustered reproductive endpoint. Do not substitute Giavi 2016 raw data or existing meta-analysis effects.
