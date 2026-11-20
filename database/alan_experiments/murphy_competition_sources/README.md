# Murphy2026 competition: unreplicated lighting evidence

Primary: https://doi.org/10.1016/j.biocon.2025.111586
Original data: https://doi.org/10.5281/zenodo.17479452
Author PDF: https://smurphlab.weebly.com/uploads/1/2/6/5/126523979/murphy_et_al._2026.pdf

Methods and Results read; the workbook and README were downloaded and their
MD5 hashes checked against Zenodo metadata. The indexed AGR record was matched
to the final DOI by title, ordered authors, journal, volume and year.

This is an eligible greenhouse setting, but lighting is **one room per arm**.
Six lamps jointly illuminate one room; pots rotate only within their assigned
arm. Thus neither pots nor lamps replicate the light treatment. The study is
retained as contextual evidence and excluded from replicated-light pooling.
No arm-level sampling variance or main-catalogue effect is calculated.

The response to competition must also be distinguished from direct lighting
responses. The original analysis reported an ALAN-by-competitor-growth
interaction for native C:N, but not for native composite growth. Competitor
growth is a measured post-treatment variable, not simply the randomized
planting density. The experiment has no zero-competitor condition. Leaf traits
do not measure herbivory or pollination, and C:N is not plant abundance.

## Original-data audit

Run `python database/alan_experiments/murphy_competition_sources/audit_sources.py`.
The script expects the external cache at
`work/gbif_annual/alan_experiment_sources/murphy2026/`.

- Plant data contain 240 unique pot keys, with ten pots in each of 24 cells.
- Chemistry has 397 rows covering 379 distinct pot-status keys after normalizing
  letter case in the documented replicate ID. All map to plant-table pots.
- Eighteen pot-status keys occur twice (36 rows); several have different values.
  Their meaning is unresolved, so no averaging or deletion is performed.
- Relative to native and exotic status for every pot, 101 keys are unrepresented.
  Missing chemistry is not mortality or absence.

`source_manifest.json` records verified source files; `pot_cell_counts.csv`,
`missingness.csv`, `chemistry_duplicate_keys.csv` and `audit.json` preserve the
audit. No full original workbook is duplicated in Git.
See the [independent review](../search/murphy2026_competition_design_review.md).
