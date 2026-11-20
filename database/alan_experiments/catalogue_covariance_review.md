# Common covariance normalization: independent check

2026-11-20. A separate agent read the ten source tables and normalized output without calling the assembler. This is a computational check, not human double extraction.

- All 108 effects and 2,006 ordered covariance entries are represented, without duplicated ID pairs.
- Values match source entries within CSV floating-point precision (maximum absolute difference 2.22e-16).
- All block diagonals match catalogue sampling variances; matrices are finite and symmetric. Minimum eigenvalue approximately -1.96e-15 is within numerical tolerance.
- Seven Lu and three Guo endpoint labels map correctly to global effect IDs.
- The Bennie/Anic shared programme retains 48 directed cross-paper entries, including 16 nonzero entries.
- The 9,658 cross-block ordered pairs are absent, not zero-filled.
- Existing effect catalogue, computed effects and effect counts were compared with Git HEAD and are unchanged.

These checks establish fidelity of normalization only. They do not resolve sampling-design limitations, publication overlap between blocks, final eligibility, or the missing study/experiment schema. No new pooled effect is estimated.

Rebuild: `python database/alan_experiments/assemble_catalogue_covariance.py`. Input hashes and block diagnostics: [audit](catalogue_covariance_audit.json).
