# Cieraad2023 primary-source design review

Reviewed 2026-11-20. [Full primary article](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.14336), Methods 2.1–2.4 and Results 3.1–3.3.

## Primary design and findings

In 2019, potted Quercus robur and Vaccinium myrtillus occupied eight treatment rows: dark, white, green and red in each of two forest types. Two posts per row supplied plant subsamples. Treatments were randomized to rows. Herbivory and leaf traits were measured after six months. Oak growth used 77 days rather than the final season because deer removed some tops; blueberry growth used 177 days. Brassica napus was lost to severe herbivory across treatments. Plants were watered as needed. The authors omitted post random effects after convergence problems.

White, green and dark treatments had similar insect damage; elevated red-light damage depended on forest type. Trait associations differed between species. Gall occurrence and abundance depended on light and forest context. No direct insect behaviour, pollination or predator rates were measured. The paper links the infrastructure to Light on Nature and the moth study [van Grunsven2020](https://doi.org/10.1016/j.cub.2020.04.083); matching individual rows remains necessary.

## Extraction implications (reviewer's inference)

- The randomization unit is the row, not a pot, leaf or post: two rows per arm overall and only one per arm within each forest type. A forest-specific effect cannot acquire sampling precision by treating the two posts or many plants as independent illumination replicates.
- Aggregate posts within rows, preserving missingness and shared dark rows. Any overall paired estimate has only two forest blocks and extremely weak variance support. Forest-type interaction estimates are descriptive without replicated treatment rows within forest types.
- Final herbivory is cumulative leaf-area loss, distinct from abundance, instantaneous consumption or pollination. Keep oak mid-season growth separate from final-time endpoints. Do not interpret associations between measured traits and herbivory as randomized mediation effects.
- Treat Dutch Light on Nature as a potential shared programme. The current article is not an independent replication of every moth report using this infrastructure.

## Original data and access

[Original dataset DOI](https://doi.org/10.34894/VD98DG), DataverseNL dataset 320747, version 1.0. [Institutional dataset description](https://pure.knaw.nl/portal/nl/datasets/replication-data-for-artificial-light-at-night-affects-plant-herb/) locates it at Klaterweg, Veluwe-Noord, near Nunspeet.

The landing page presents an Anubis challenge, but the public [metadata API](https://dataverse.nl/api/datasets/:persistentId/?persistentId=doi:10.34894/VD98DG) returned JSON successfully. Raw workbook download and checks are handled separately. Do not record the dataset as inaccessible solely from its landing page.

The publisher lists `jpe14336-sup-0001-AppendixS1.pdf`; [supplement URL](https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F1365-2664.14336&file=jpe14336-sup-0001-AppendixS1.pdf) returned 403 in this review. Exact spectra/intensity and post coordinates remain to be verified from the supplement or raw metadata; do not invent them.
