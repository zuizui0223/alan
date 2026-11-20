# Independent methods audit — 2026-11-20

Scope: read METHODS_ja.md, the three synthesize_*.py implementations, synthesis READMEs, Cieraad/beech extraction READMEs and search/README.md. No new analysis or literature retrieval; METHODS_ja.md was not edited.

## Verdict

No numerical-method contradiction found. Capture synthesis correctly averages raw cell means before logs (eight illuminated cells and two dark cells in SevenPonds); aphid and leaf syntheses average log contrasts. Counts of programmes are 2/2/3, with 4 dependent aphid and 11 dependent leaf contrasts. Covariance and a'Va aggregation, k=2 analytic REML checks, modified Knapp–Hartung max(1,q), and df=k-1 match code. Penryn first-report final observation is distinguished from programme-terminal follow-up. The draft correctly avoids completed-systematic-review, preregistration, human-double-extraction and mechanistic-mediation claims. Existing habitat intervals remain a separate exploratory within-pond comparison.

## Precision improvements before publication

1. Leaf-damage paragraph: add “retained sampled leaves/plants” and identify damage as fractions/percentages, not absolute consumed area or whole-plant loss. Source READMEs and script explicitly retain attrition uncertainty; Swiss all-leaf mean percent is not total consumed area divided by total leaf area. Current “葉面積被害量” is unnecessarily ambiguous.
2. Aphid paragraph: include the measurement-compatibility assumption stated in its README: Penryn uses three-minute search counts whereas Miller uses cage counts; interpreting their ratios together assumes treatment-relative detection differences do not generate the contrast. The current paragraph states neither observation distinction nor that assumption.
3. Search paragraph: clarify that the ALAN term group was restricted to title/abstract, rather than implying every search term was title/abstract-restricted. The search README describes this specifically. Retain the incomplete WoS/Scopus coverage and evolving adjudication caveats.

These are bounded reporting clarifications, not reasons to recompute or broaden the synthesis. The broad statement that selection was independent of direction/significance is a workflow assertion, not something the three numerical scripts alone verify; it should remain linked to the selection/plan records.

Main-agent resolution: all three wording clarifications incorporated into METHODS_ja.md on2026-11-20. No effect values or selection changed.
