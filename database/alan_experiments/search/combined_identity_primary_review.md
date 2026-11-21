# Primary bibliography adjudication of three supplementary-search candidates

2026-11-21. Metadata were retrieved from the official EuropePMC core-search
API and bioRxiv API. Source response hashes and retained bibliographic fields
are in `combined_identity_primary_metadata.json`. This resolves identity,
not experimental eligibility, full-text reading or numerical equivalence.

## Same article, two index records

- [AGR IND601326109](https://europepmc.org/article/AGR/IND601326109) and
  [MED25914438](https://europepmc.org/article/MED/25914438): title and ordered
  author surnames match; Ecological Entomology40(3),187–198, June2015 match.
  MED supplies DOI10.1111/een.12174, absent in AGR. Online first-publication
  dates differ from the print issue date; this is not evidence of another
  study. Initial abbreviations vary, but do not contradict ordered authors.
- [AGR IND605878087](https://europepmc.org/article/AGR/IND605878087) and
  [MED27251575](https://europepmc.org/article/MED/27251575): matching title,
  ordered authors, Global Change Biology23(2),697–707, February2017.
  MED supplies DOI10.1111/gcb.13371. Online and print dates again differ.

The article landing pages failed in the web reader, but the official
`https://www.ebi.ac.uk/europepmc/webservices/rest/search` endpoint returned
the four complete core metadata records using `EXT_ID:ID AND SRC:AGR/MED`,
`resultType=core`, `format=json`. No full-text reading is inferred.
The two resolutions are explicit in `combined_search_identity_resolutions.csv`;
the original search snapshots and source DOI blanks are not overwritten.

## Preprint and journal version, not identical data

[Official bioRxiv metadata](https://api.biorxiv.org/details/biorxiv/10.1101/2023.07.03.547605)
for version1 explicitly sets `published` to10.1038/s41598-024-73378-9.
The title and ordered five authors agree with the
[final article](https://www.nature.com/articles/s41598-024-73378-9).
This verifies a publication-family relationship. It does not prove identical
methods or observations: the final article includes later experimental
conditions. Neither report is converted into an additional independent
experimental replication or assigned an effect here. An independent AI
reader confirmed the official version link and author correspondence.

This additional link is stored separately from the original-search-only
version table, whose domain excludes some supplementary records.
`combine_search_snapshots.py` joins both verified-link tables for the combined
inventory. The result is546unique record IDs,539reports and535publication
families. The earlier541report-key count preceded these two identity
resolutions;572search appearances and the26record-ID overlap are unchanged.
The three candidate pairs are resolved, but this does not assert that all
possible duplicates across all searches or all experimental programmes have
been resolved. Final selection and the planned subscription-index searches
remain incomplete.
