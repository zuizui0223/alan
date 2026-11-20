# Ceppi2025: beech herbivory and phytochemistry

Reviewed 2026-11-20. Verified DOI [10.1111/1365-2664.70040](https://doi.org/10.1111/1365-2664.70040). [Publisher Methods 2.1–2.5 and Results 3.1–3.4](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.70040) read in full.

## Primary design and findings

At two Swiss forest sites,56 potted beech trees received dark,21,35 or70lux treatments for43 nights in2022. Seven trees per site/treatment stood in a quarter-circle around their light source. Leaf morphology, chemistry and chewing damage were assessed at the end. Damage probability decreased only under70lux; damage severity among already-damaged leaves did not vary by lighting. Some phytochemical responses were non-monotonic. The phytochemistry–damage associations were correlational. No predators or pollination were measured. Sites were Birmensdorf and Alpthal, and funding identifies ALANeX. Original data are [EnviDat10.16904/envidat.628](https://doi.org/10.16904/envidat.628).

## Extraction implications and original-data access

Seven trees and15 leaves/tree are subsamples of a shared light installation, not independent lighting replicates. The described design implies one installation per dose per site; verify installation identifiers/layout before using precision. At most two site blocks support treatment contrasts from this design, and dose arms share their controls. Damage incidence and severity conditional on being damaged are distinct endpoints; neither can silently replace unconditional leaf-area loss. Treatment-associated chemistry is not a randomized mediator intervention.

The Swiss sites differ from Cieraad's Dutch experiment. Audit potential ALANeX infrastructure overlap with van Koppenhagen2024 before claiming independence from other Swiss reports; shared funding alone does not prove shared plots.

The [public metadata API](https://www.envidat.ch/api/3/action/package_show?id=raw-data---artificial) returned the downloadable workbook, which was retrieved successfully (197472bytes; SHA256 `b7c0c7a81315db06951cabcfcb070f926236a42d1f05c0b93b65333cad88e78d`). Metadata and workbook are cached externally under `alan_experiment_sources/beech2025/`. Chemistry sheet contains56 rows; leaf-damage sheet855 rows for57 tree IDs, each with15 leaves. The extra tree is in Alpthal treatment `2200_100`; the stated experiment's four groups each have seven trees/site. Preserve this extra treatment separately rather than pooling it with the nominal100% arm. These are source row counts, not independent lighting sample sizes. No effect or variance has been extracted in this review.
