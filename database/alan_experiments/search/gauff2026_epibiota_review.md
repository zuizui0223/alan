# Gauff2026 seagrass–epibiota review

Reviewed 2026-11-20. Read Methods/Results in [institutional pre-proof PDF](https://archimer.ifremer.fr/doc/01074/118553/132576.pdf), DOI [10.1016/j.marenvres.2026.108273](https://doi.org/10.1016/j.marenvres.2026.108273), and all tables in [supplement](https://archimer.ifremer.fr/doc/01074/118553/132577.docx). Direct requests returned HTTP 200, respectively 1206454 and 39221 bytes; the web tool PDF request returned 403. No raw observation table inspected here.

## Primary-source summary (under 200 words)

Six alternating Venice Lagoon plots received ALAN or control conditions, three per arm. Each contained six quadrats; three were cleaned of epibiota every two weeks. Thus removal was repeated in the field, not merely an immediate assay manipulation. Exact collection timing relative to each cleaning was not specified. Night-time measurements were made in situ; laboratory shoots underwent at least eight hours of darkness, then fifteen minutes of standardized daylight.

Night-time oxygen flux and fluorescence showed no detectable treatment effects. Laboratory daylight quantum yield and leaf length showed interactions; biomass did not. However, the leaf-length post-hoc comparison is inconsistent: text attributes significance to cleaned lit versus cleaned dark, whereas supplementary Table 11 gives p=.2574 for that contrast; p=.0027 compares uncleaned versus cleaned within lit plots. Supplementary oxygen Tables 2 and 3 are identical despite different endpoint labels. Tables 7/8 have inconsistent dark/light-adaptation headings. These discrepancies require resolution before numerical interpretation.

## Independent design and extraction implications

- Independent lighting replication is three plots per arm, with alternating placement rather than verified randomized assignment. Quadrats, sampled shoots, chamber replicates and repeated dates cannot increase lighting n. Quadrat-level epibiota manipulation is nested within the lighting plot. Preserve within-plot covariance for treatment modification.
- The removal treatment is documented as biweekly from permanent quadrats. Do not call it only instantaneous cleaning. Conversely, the paper does not identify the exact delay between the latest removal and each collected shoot or isolate physiological effects of mechanical cleaning. A whole-treatment contrast represents repeated epibiota removal with its associated handling, not a proven epibiotic mechanism alone.
- Standardized laboratory daylight tests are responses after field exposure history; they are neither independent lighting experiments nor direct measurements of in situ daytime carbon gain. The dark-adaptation and daylight values share exposure plots and sampling history.
- Field oxygen outcomes can be negative. Ordinary lnRR is not defined for nonpositive means and must not be rescued by arbitrary constants or absolute values. Production is a difference of oxygen measurements from two different shoots in the same quadrat, not two measurements on the same shoot; its uncertainty needs correct pairing and measurement provenance.
- Fluorescence efficiencies, leaf length, shoot biomass and SLA are distinct endpoints. A significant interaction in one cannot establish increased net production, population growth, predation or pollination. A three-plot arm mean does not gain precision from the 33–65 observations listed in model tables.

## Independently verified supplementary discrepancies

| Table/comparison | Estimate | SE | Statistic | Reported p |
|---|---:|---:|---:|---:|
| S11 cleaned lit − cleaned dark | −24.084 | 13.311 | −1.809 | .2574 |
| S11 uncleaned lit − cleaned lit | 31.981 | 9.233 | 3.464 | .0027 |
| S8 cleaned lit − cleaned dark | 0.0969 | 0.0356 | 2.7163 | .0328 |
| S8 uncleaned lit − cleaned lit | −0.0944 | 0.0282 | −3.3424 | .0046 |

S11's two contrasts are different scientific questions. Do not relabel the significant within-lit cleaning contrast as an ALAN effect among cleaned quadrats. S10 still reports a leaf-length interaction, so this mismatch is not grounds for erasing all reported interaction evidence. S2 (net production) and S3 (gross production) repeat all displayed fixed estimates, random-effect entries, observation count and R². Their numerical identity is verified, but the correct table cannot be inferred. S7's title specifies light adaptation while its outcome header says dark; S8's caption says dark despite the article using it for the daylight assay. Preserve source labels and flag the ambiguity rather than silently correcting them.

## Decision

Field factorial experiment is eligible in principle for a seagrass physiology/epibiota-modification stratum. Keep reported model contrasts separate and on hold pending raw data, explicit coding, quadrat/date mapping and table reconciliation. No self-extracted effect size, reconstructed arm SD or causal mediation estimate was produced. This is marine epibiota manipulation, not direct evidence for freshwater wetland plant decline or a terrestrial predator cascade.
