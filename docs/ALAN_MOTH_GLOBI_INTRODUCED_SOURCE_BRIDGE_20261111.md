# ALAN × moth-like flowers × island naturalization: real-source feasibility, 2026-11-11

**Executed real-source audit, not a fitted ALAN or predation effect.** This document accompanies draft PR #40 and the reproducible signed-artifact source join. Neither fd nor island Chapter 1 was edited.

## Executed sources and commands

Read-only authorized ZIPs:
1. island Actions run **34232450884**, artifact **10058653212**, containing September 2026 All/Direct trait snapshots and historical island-specific origin ledgers.
2. mois Actions run **41976552977**, artifact **13391803415**, containing **61,967** GloBI positive butterfly/moth flower links from the bounded November 2026 source.
3. mois Actions run **42171173575**, artifact **13481850056**, containing 889-island ROI-free original photograph observation support (2012–2024).

Run with actual downloaded artifacts, without exporting protected species rows:

    python scripts/analysis/alan_moth_source_bridge.py \
      --island-artifact /secure/island_artifact.zip \
      --globi-artifact /secure/globi_stage2.zip \
      --photo-artifact /secure/photo_source.zip \
      --output /secure/aggregate_result.json

The script asserts five source-member SHA256 hashes, 61,967 authentic interaction rows, unique island×species status and photo island×year keys, and valid categorical sources. A full local real-data execution completed; see results/alan_moth_source_bridge_20261111_public_summary.json. Existing outputs and scripts were not overwritten.

## Empirical support and the decision it forces

- 15,382 species names have simultaneously recorded colour and flower form (All scope), and 1,918 have all three of colour, form, and tube depth.
- EXACT white + one tubular/salverform/funnel_trumpet state gives **950** candidate names; **881/950** lack tube depth. Of 69 depth-known candidates, **10** have documented deep tubes.
- The ten deep-tube candidates have **zero matched GloBI moth flower-visit records**. This is absence of database documentation, not evidence of biological absence or failed pollination.
- Among the broad 950 candidates, 14 have a positive GloBI moth visit (4 moth-recorded-only, 10 both GloBI moth and butterfly); no study demonstrated exclusive nocturnal dependence.
- Historical island-specific flora includes 7,477 candidate island×name rows on **1,740** islands; only **183 rows / 48 names / 67 islands** carry an *introduced* origin code. The remaining 6,027 candidate island×name rows have unresolved origin; introduced does not certify established/naturalized.
- Crossing broad shape + positive GloBI moth record + the historical introduced origin gives only **5 distinct names / 8 islands** (10 island×name source rows), BEFORE adding light or geckos. These are not gecko/light-complete island samples.
- Historical Direct-only sensitivity: **738** broad candidates; **683** lack tube depth; **10** strict deep; **67** islands with introduced broad candidates and **8** with additionally positive GloBI moth records.

The November rebuilt photo-island polygons and September historical status IDs appear compatible in text, but textual IDs are NOT sufficient proof of same GSHHG physical geometry: counts involving 889 photo-island ID intersection are **provisional only**, not a verified spatial join.

No light × gecko × naturalized flower three-way interaction is estimable from the linked source panel. In particular, 8 sites is not a viable independent archipelago population, and zeros for undocumented statuses may not be inserted.

## Better naturalized flora source: GloNAF 2.02, not GRIIS country checklist

Public GloNAF v2.02 (Davis et al. 2025, Ecology, DOI 10.1002/ecy.70245, openly hosted at https://zenodo.org/records/17105725) lists **16,429 harmonized naturalized/invasive taxa and 1,343 regions, including 427 island regions**. It is the correct independent input for a genuinely naturalized endpoint.

Acquisition contract for the next executed join:
1. Pin Zenodo record **17105725** and published file MD5s. Download glonaf_flora2.csv, glonaf_taxon_wcvp.csv, glonaf_region.csv, glonaf_list.csv plus glonaf_2024_regions.zip; record bytes/hashes and CC BY 4.0 attribution.
2. Join glonaf_flora2.taxon_wcvp_id to glonaf_taxon_wcvp.id and use WCVP accepted_plant_name_id where complete; never count the entire taxon table independently of active flora rows. Preserve original name and taxonomic match quality.
3. Require glonaf_region.island == 1 AND inspect the actual polygon. Do NOT assign a multi-island region list to every GSHHG member island. Keep composite archipelago records at their legitimate coarser spatial unit or explicitly exclude them from single-island modeling.
4. Keep source-list coverage and completeness fields, source date and taxon disambiguation. Unlisted candidate species are **NOT established absences** unless the underlying checklist is assessed as near-complete for naturalized vascular plants.
5. Re-estimate observed-and-expected genus/source-conditioned proportions on the *naturalized* GloNAF cohort, NOT on historical island introduced status or raw GIFT native tables.

At the time of this PR, the GloNAF actual files are **not yet downloaded/joined with the archived private plant traits**; **427 is candidate published region coverage, not 427 analysis-ready geocoded island units**.

GRIIS Country Compendium explicitly notes that country-level records normally pertain to mainland units and cannot confirm offshore island gecko nonpresence. An island-specific reptile census/range source is needed. GBIF or iNaturalist positive gecko detections may supplement **presence only**, never confirmed absence.

## Annual light and temporal inference

EOG VIIRS VNL versioned ANNUAL grids (2012–2024, 15 arc-second, radiance plus cf_cvg) are the target source; see https://eogdata.mines.edu/products/vnl/ and Earth Engine catalog NOAA/VIIRS/DNB/ANNUAL_V22. A 2025 image should not be labelled V22 without a genuine source product. Spatially intersect verified polygon geometry, keep cf_cvg and satellite/product/year, and use stable area/within-island covariates. Do not substitute 2020/2024 quarterly proxies for annual values.

Current naturalized checklists and current nightlight are cross-sectional: plant establishment may predate artificial lighting. Even a significant association cannot prove that moth-dependent plants were prevented from establishing. A temporal estimand requires credible documented establishment dates, or time-matched photo/occurrence observations with explicit observer effort.

## Analysis choice

**Promote the source-rich observational arm** to primary feasibility exploration:
- Spatial units: island × repeated site × observation year (2012–2024), only real original photo UUIDs and accepted species.
- Trait: broad *white + tubular/salverform/funnel_trumpet* as a fixed, explicitly proxy phenotype; "deep" only when independently observed, unknown kept unknown.
- Endpoint: change in *photographed taxon-type share*, with equal accepted-species weighting, observation/photographer opportunity, detection, urbanization, source composition, season, and effort controls. This is **not a flora census or plant decline**.
- Cross-sectional naturalized GloNAF with complete, independently reviewed gecko statuses is a separate secondary analysis.
- Mechanistic validation: independent pollen deposition/seed set plus concurrent nightly moth and gecko behavioral observations; neither GloBI moth visit nor VIIRS upward radiance measures that chain.

Avoid optimizing for significance. Minimum viable support is real geographic overlap and independent archipelagos with verified gecko presence **and verified absence**; conduct spatial block holdout and simulation-based power analyses before any 3-way coefficient. Do not treat 8 positive island records, 889 photo islands, or 427 GloNAF regions as independent replicated contrasts for the triple-interaction.

## Validation and privacy

Synthetic tests: 7/7 passed locally, including source bridge cases. No protected species-to-trait or species-to-island rows from the 100 MB private artifact are committed. Published result is aggregated; original source files remain with their owner/source providers. Claims and historical results in existing ALAN, fd and island documents are unchanged.
