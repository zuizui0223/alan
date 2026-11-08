# ALAN — 512 new independent island photo observations to densify site-year support (2026-11-08)

**Executed source-only outcome: [GitHub Actions #42454984788](https://github.com/zuizui0223/alan/actions/runs/42454984788), SUCCESS.** The source target file contains **512 genuine original observation UUIDs and exactly one original photo ID per UUID**, all unmeasured in the first 1,280 source photos (batches 1–5). It is not a set of newly downloaded images or a claim about flower colour.

## The actual source-bottleneck analysis

The 8,942 SHA-locked new island photograph opportunities span 4,776 genuine GSHHG island × 0.05° site × original observed year cells. But the current island-breadth acquisition order initially spreads photographs over many singletons and cannot immediately create many cells that are eligible for the existing 5-photo/3-WCVP-accepted-species ALAN strict estimator.

Holding the original 1,280-source-photo batches (1–5) fixed, the same **additional 512 authentic independent observation UUIDs** were considered under two predetermined pre-colour routing rules:

| Pre-colour source routing with 1,280 first-batch photos held fixed | Number of original source site-years with at least 5 original observations and 3 observed accepted taxa |
| --- | ---: |
| Only first 1,280 source photo opportunities | 3 |
| Next 512, breadth-first island schedule (1,792 photo opportunities total) | **6** |
| Next 512, **site × observation-year densification** (also 1,792 photo opportunities total) | **60** |

**This tenfold comparison concerns source PHOTO OPPORTUNITY ONLY, not the actual number of colour-classifiable strict ROI9 sites.** Even a source-qualified 5-photo/3-taxon cell is not a valid *colour-measured* cell if its photographic ROIs fail. It does not mean that 60 new strictly eligible ecological sampling cells already exist.

The targeted 512 original source photo IDs are distributed across **20 independent GSHHG islands and 61 original site × year cells**. The route is deliberately different from the island-breadth source: island-wide source support is initially weaker than geographically broad source by design, but photographic trials are consolidated at places where actual year-specific VIIRS and WorldPop values were already verified.

## How the source was chosen

1. Read the same immutable source-only 8,942 observation/photo schedule (SHA256 `e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd`) and require 8,942 source observation UUIDs with unique source photo IDs. The source record itself contains **no new photo-colour measurements**.
2. Remove all 1,280 source photos/observations in original batches 1–5 (including the 1,024 images processing under [PR #7](https://github.com/zuizui0223/alan/pull/7)). Source UUIDs cannot be placed in a new colour worker twice even if prior photo inference is still running.
3. Find original GSHHG site × observation-year cells with at least **five distinct genuine original source observations and three distinct accepted plant taxa** in the complete 8,942 candidate frame. Exactly **301** uncompleted sites had qualifying source opportunity.
4. Prioritize islands with zero previous strictly source-qualified ALAN photographic cells, then potential source support and accepted-taxon breadth, then SHA-based colour-blind source ID ordering. Limit at most **64 photos per island and 16 per site × year**, selecting distinct accepted taxa before different observations of the same taxon.
5. Freeze **512 true new observations**, produce two disjoint 256-photo source work packets, and maintain original dates, coordinates, source taxon identifiers and year-specific environmental availability. No image segmentation or colour inference occurs at this step.

## Statistical and ecological interpretation

The original global source photo index is not a census of confirmed wild-native flowering plants. The new 512 are adaptively chosen by known landmass/source support, so **their source sampling probabilities are not the original 93m tier06 random-hash probability**. Do not combine the 512 with the historical-main or island-breadth sample to estimate an unweighted worldwide or island-wide floral whiteness rate. In a later regression the data need explicit source-design stratification and observation/taxon/site clustering.

The photographic-ROI9 colour success rate in the independently completed first256 island-breadth run was 92/256, not 256/256. A larger source-opportunity-qualified denominator does **not** remove the actual need for three distinct successfully measured flower-region ROI species. This new route is a way to raise *the number of opportunities* to achieve that condition, and its effectiveness must be measured using actual original ROI outcomes.

## Signed files

[Source-only Actions #42454984788](https://github.com/zuizui0223/alan/actions/runs/42454984788) is SUCCESS. The ephemeral Actions artifact `alan-original-observation-512-site-densification-source-no-colour` ID **13606624937** contains `alan_site_year_densification512_original_photo_source.csv.gz` SHA256 `172680e798ba9711ee149fd75fc84601e55137c17a818c028aa36ea24c309697`, its per-batch source-coverage aggregate, and SHA source receipt. Numerical source-only summary is also in [results JSON](../results/alan_site_densification512_original_source_20261108.json). The source selection code and tests are retained under `scripts/acquisition/alan_site_densification_from_8942_20261108.py` and `tests/test_alan_site_densification_from_8942_20261108.py`.

**MOIS remains moisture-only.** This ALAN source-only extension does not alter any pre-existing five-year island×nightlight model coefficients or independent original-photo colour measurements.
