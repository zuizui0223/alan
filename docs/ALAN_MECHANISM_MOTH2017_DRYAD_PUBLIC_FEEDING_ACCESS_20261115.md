# ALAN mechanistic light-sensitivity source: original moth feeding Dryad dataset availability (2026-11-15)

## Original candidate

van Langevelde, van Grunsven, Veenendaal & Fijen (2017), Biology Letters 13:20160874, DOI 10.1098/rsbl.2016.0874; original open Dryad record [Data from Artificial night lighting inhibits feeding in moths](https://doi.org/10.5061/dryad.8v212).

The authors tested **four moth species** under dark, red, green and white light with individual moth IDs, experiment dates, blocks and observations of feeding activity every roughly **six minutes**. The source field labels advertised by Dryad are Day, Block, Treatment, Species, Sex, unique.id, Obs and Size. A per-observation feeding indicator is not an independent moth, an effective pollen transfer event, or a reproductive plant response.

### Actual public source acquisition result

Executed [ALAN GitHub Actions #43706153523 SUCCESS](https://github.com/zuizui0223/alan/actions/runs/43706153523). Success = four parser/unit tests passed and a reproducible endpoint receipt was generated, **NOT download of original individual measurements**.

- Dryad anonymous official API dataset DOI metadata HTTP 200, versions metadata 200; one version files metadata HTTP 200.
- Original public web file stream path 27006 returned HTTP **403**.
- Dryad API file download paths returned HTTP **401** in this environment.
- No individual source moth rows were obtained, no observations analyzed, **no numerical per-species or per-spectrum source effect is claimed**.
- Repeated public-source metadata links may refer to other dataset identifiers. Such links are not new moth experimental studies and are not pooled.

Reproducibility: scripts/analysis/alan_mothfeeding2017_Dryad_source_20261115.py, tests/test_alan_mothfeeding2017_Dryad_source_20261115.py, receipt artifact #14157788426 (ZIP SHA256 77c3663d4208241a420f738ede6fba02f3095d5cac90f4c3684239499958905a). The script refuses to fit source effects until original source column schema, individual IDs and matched light labels are verified.

## Correct analytical role

The paper supplies external evidence that moth **feeding behavior** can depend on lighting and spectrum, and—if its individual data become retrievable—could help quantify a *behavioral nocturnal guild light-sensitivity* modifier. But it does NOT measure effective pollen arrival, pollinator visits per flower with verified pollen transfer, fruit set or plant seed production, and cannot by itself provide the central A slope on night dependence. Do not infer the loss of plant reproduction as feeding decline × inferred moth syndrome.

The original ALAN review's strict candidate count remains 12 plant reproductive papers, 9 direct ALAN manipulated; this 2017 moth behavioral paper is a separate mechanistic source, never an independent plant-reproduction experiment.
