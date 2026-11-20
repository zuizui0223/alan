# Giavi2020: reproduction and pollinating seed predation

[Primary paper](https://doi.org/10.1038/s41598-020-68667-y), Methods/Results/Table1 and public supplement checked2026-11-20. This report was already primary-reviewed in batch1; the update adds extraction, not a new reviewed paper. The supplement checksum matches the primary XML; it contains Matlab video-processing code, not field observations or fitted model objects. It was inspected as text and not executed.

Seven field-model coefficients with SEs are retained in `reported_field_coefficients.csv`. Six are transcribed from Table1; the illuminated-versus-dark infection coefficient is explicitly reported in Results. It agrees with the difference of the other two infection coefficients, but its SE is taken from the reported comparison, not constructed by assuming zero covariance. All seven remain outside the group-mean lnRR catalogue.

The seed-development outcome for Epilobium is whether a selected capsule produced at least four seeds; it is not whole-plant fruit set or pollen deposition. Silene selected-capsule seed development and seed numbers in both species have reported nondetections without numeric treatment coefficients in this table; retain those outcomes as unextracted, not zero effects. Infection refers to Hadena bicruris seed-predator larvae in capsules, not adult-moth predation. Hadena's pollinating role does not make this a direct estimate of its per-visit pollination efficiency.

## Model interpretation holds

Eight meadows comprise four lit and four dark controls. Adjacent dark sectors share the four illuminated meadows. The paper also calls these sectors sites; whether the fitted site random effect represented eight meadows or12sectors cannot be certified without code. Plant and capsule observations are not lighting replicates. Relationship to the other Swiss meadow programmes remains unresolved.

The stated relative-reproduction endpoint is intact capsules divided by exposed flowers, but the binomial Methods describe two response columns as intact capsules and total exposed flowers. If implemented literally, those columns imply intact/(intact+flowers), not the stated proportion. This may be a writing error or an implementation issue; neither is assumed. Preserve the original coefficients with a denominator hold. Do not derive a risk ratio, a lighting-versus-dark reproduction SE, or mediation from them.

The link is not explicit in the paper and is not silently assigned from a software default. Estimates remain on the author-model scale. No model refit or pooled reproductive estimate is claimed. Three infection contrasts are dependent recodings of one outcome, not three replications. A separate reader independently verified all seven coefficients and the definition/dependence issues.

Reproduce: `python database/alan_experiments/giavi2020_sources/extract_reported_coefficients.py`. Source hashes and acquisition provenance are in `source_audit.json`. Public field data remain author-request; no author contact was made.
