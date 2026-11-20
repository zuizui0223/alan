# Bennie2018: Lizard grassland, primary data recovered

Reviewed 2026-11-20. Full Methods and Results read from the [institutional primary PDF](https://researchmgt.monash.edu/ws/portalfiles/portal/248033872/65639617_oa.pdf), DOI [10.1111/1365-2664.12927](https://doi.org/10.1111/1365-2664.12927).

At Lizard National Nature Reserve,30 plots received five treatments, six plots each: high/low white, part-night low white, amber and unlit. Lighting began April2012. Cover was measured by18 point-quadrat pins annually2011–2015; species-sorted dry biomass was harvested annually2012–2015. Fortnightly flowering observations generated inflorescence-weighted mean dates. The field had poorly drained gley soil and twice-yearly mowing. Total biomass and annual flower density showed no significant treatment effects; species-level biomass/cover and flowering responses differed among years. This natural grassland differs from the54 Penryn artificial mesocosms, although it used related lighting technology.

## Data access resolved

The paper's [Exeter handle10871/27138](https://hdl.handle.net/10871/27138) redirects to Figshare article29731352. Its landing page returns202/empty, but the [official metadata API](https://api.figshare.com/v2/articles/29731352) returns the archive. [File56750435](https://ndownloader.figshare.com/files/56750435) was downloaded:59,084bytes, verified MD5`6332d8b00cf937916248dc68c7cdd230`. It contains15files:12species-year phenology CSVs, `biomass_v3.xlsx`, `harvested_weights.csv`, `point_quadrat_data.xls`. Cache:`alan_experiment_sources/grassland2018/`.

## Extraction inference

Recommend **eligible_field_design_original_data_acquired_pending_extraction**. Independent units are plots; pins, species and repeated years are not new experiments. Preserve shared dark controls and repeated-plot covariance. Final-year biomass can be extracted separately from cover and phenology after checking IDs and units. Weighted mean flowering date is not flowering peak, flowering duration or animal–plant mismatch. Use day differences, not lnRR of calendar dates. Annual sums of observed flower counts are not counts of unique flowers or seed production. Do not call this poorly drained grassland a dry-system control. Exact overlap with other Lizard papers still requires a plot crosswalk; do not combine its programme with Penryn merely because authors and hardware overlap.

No numerical effect extracted in this review; original ZIP files remain unmodified.

## Workbook schema audit (original data, no effect calculation)

`biomass_v3.xlsx` has one sheet,120 rows and19 columns: `Plot`, `Light`, `Year`,15 material/species columns including `Dead `, and `total`. There are30 unique plots, each once in2012–2015, six per arm-year, no missing cells, and invariant treatment per plot. Codes are `A`, `CON`, `HIW`, `LIW`, `LIWT`. The final biomass year is2015; exact harvest dates are absent, while Methods specify September harvesting. Final-year totals match the sum including dead material. One earlier row, G7/HIW/2012, differs by0.125 from the component sum.

`harvested_weights.csv` has29 rows and no year field. On its29 matched plots, the three dominant grass columns and total equal exactly four times their2012 workbook values. Unit/scaling documentation must therefore be checked before reporting absolute biomass; do not append this as an additional harvest or silently replace the workbook.

Phenology schemas vary: later files use `stage1_3`, `stage4_6`, `stage7_9`; Holcus2012 instead has nine stage columns, some blanks and `#REF!` density entries. Recorded sample areas vary within2012 (1,0.25,0.06) and are absent from Holcus2012, whereas inspected later files use0.06. Do not assume these columns already represent the1m² flower totals described in Methods. `point_quadrat_data.xls` needs an XLS reader; it was not parsed here because `xlrd` is absent.

The external cache also contains `archive_inventory.csv`: all15 ZIP member paths, byte sizes, CRC32 and SHA256. No archive member was modified or extracted onto a repository path.
