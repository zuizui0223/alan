# EERES arrangement: verified original workbook and compensation decomposition (2026-11-25)

User supplied the public Dryad archive DOI 10.5061/dryad.rfj6q57k5, version 20251009. Archive SHA256 2f16f4aa160847aa23132a0be73dde8e73131cd7510d102675dd631b9440f769. Inner DISCONNECTED_alldata.xlsx SHA256 496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f, exactly agrees with prior source lock. README confirms 18 network flumes in six shared header blocks, 26-day colonization and 23-day treatment. Original workbook read locally; 144 drift records/5564 individuals (38 taxa), 108 benthos records/11326 individuals (40 taxa), 324 background rows. Verified treatment identity, two upstream sides per flume/week, four weeks, two final downstream benthic subsamples/flume. Raw file not redistributed.

## Actual network comparison (resolves previous bounds)

CT two untreated branches, SS L and F separately, CS L+F and one untreated. Sum the two upstream trap vectors within flume-week, count individuals and richness of taxon union, average four weeks, paired across six HT. This is summed tributary capture, not a measured downstream flux. Downstream benthos measured separately at final sampling; average two subsamples per flume. Ratio of block-average means, paired delta t5 intervals. Source Python scientific computation executed in the current session; block vectors and core estimates independently recalculated in JS. Saved replay adaptation uses pandas reader and checks original SHA. No CI execution in this session.

| Endpoint | CT | SS | CS | CS/SS | approx 95% CI | nominal p | Holm four CS/SS tests |
|---|---:|---:|---:|---:|---|---:|---:|
| drift_count | 95.542 | 74.583 | 61.708 | 0.827 | 0.689–0.993 | 0.0447 | 0.1787 |
| drift_rich | 12.167 | 12.292 | 11.375 | 0.925 | 0.824–1.039 | 0.1458 | 0.4236 |
| down_benthos_count | 94.167 | 127.667 | 121.000 | 0.948 | 0.806–1.114 | 0.4336 | 0.4336 |
| down_benthos_rich | 10.333 | 10.250 | 9.250 | 0.902 | 0.776–1.050 | 0.1412 | 0.4236 |

Primary descriptive result: CS/SS drift count .827 (5/6 blocks negative, leave-one-out .781–.859); raw count totals CT2293 SS1790 CS1481 sum5564. Sign-flip reference p=.09375 (not actual randomization test). SS/CT .781 [.696,.875] and CS/CT .646 [.548,.761]. Arrangement effects on richness and downstream standing benthos not established. Direct drift-vs-benthos count ratio of RRs .873 [.729,1.045]; not significant, don't contrast significance by itself. Temporal CS/SS point ratios .952,.790,.827,.761; T4/T1 .800 [.525,1.218], all post-treatment, not a confirmed time trend.

## Additive-scale identity, not mediation

Let C mean untreated branch in CT flume, X untreated branch in CS flume, L,F,LF other branch means. The exact within-block identity CS−SS = [(LF−F)−(L−C)] + [X−C]. This partitions the observed difference into a local interaction-shaped contrast and untreated-context contrast, not a separate randomized test of context-independent local effects. In this design, local treatments occur in different flume arrangements, so background/neighbor effects may enter the first term too.

Observed C=47.771, X=48.625, L=59.583, F=15.000, LF=13.083. X/pooled-three-CT mean=1.01185, within previous 1.27977 sign threshold. X/C=1.01788 [0.8482,1.2215]; no evidence of compensating increase or equivalence. Local term −13.729 [−24.532,−2.927]; context term +.854 [−7.935,9.644]; total −12.875 [−25.334,−.416]. Arithmetic identity supports a candidate explanation: lighting-related drift increase on flowing branch disappears in ponded branch, while unperturbed branch does not show large compensation. Not evidence that actual insect export to land or pollination is reduced. No claim of equivalence of control branches or causal mediation percentage.

## Scope and provenance

Previous partial identification file retained as history; now exact branch identity resolves it. This is one programme already in PR101 and earlier stage analysis, not new independent k. Holm correction within four newly computed CS/SS endpoints only; exploration wider than this family. No definitive general effect, no macro plant causal mechanism, main untouched. The repository JSON carries the source-derived block vectors, contrasts and decomposition, and week-summary results from the local source computation. 
