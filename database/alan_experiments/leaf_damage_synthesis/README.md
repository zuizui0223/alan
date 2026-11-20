# Exploratory meta-analysis of terminal leaf damage

Reproduce: `python database/alan_experiments/synthesize_leaf_damage.py`.
Selection and estimand were recorded in [the plan](../leaf_damage_meta_plan.md)
before the pooled estimate was computed.

Three independent experimental programmes contribute 11 dependent contrasts:
Netherlands, six species-colour contrasts; Switzerland, three doses; Hong Kong,
two plant endpoints. Within-programme log contrasts are averaged equally with
their full sampling covariance. Species and lighting menus differ, so this is
not a standardized common-dose effect. No wet-versus-dry comparison is made.

The programme response ratios are 1.908, 0.680 and 0.720, respectively.
REML with modified Knapp-Hartung uncertainty yields a pooled ratio of 0.929
(approximate 95% interval 0.237–3.633; k=3, df=2; tau-squared=0.210).
The result does not establish a common increase or decrease in damage.
It also does not establish absence of an effect or equivalence.

Hong Kong is newly aggregated to four sites, with paired plots averaged
within site and lighting arm. This is different from the previous seven-pair
estimate; both are preserved, not silently substituted. The two plant
contrasts have strongly negative covariance, producing a small variance of
their mean. This cancellation over four sites is not proof of precise or
homogeneous responses across plant species. Netherlands and Switzerland each
have only two cluster pairs. These small-cluster delta variances remain
uncertain despite the small-k meta interval.

The response is damage in retained sampled leaves/plants, not attack rate,
herbivore abundance, whole-plant loss, survival, or population decline.
Leaf aggregation and observation duration differ. The Swiss endpoint includes
undamaged leaves; its conditional-severity and incidence endpoints are excluded.
No other report from either the Dutch or Swiss infrastructure is included,
avoiding duplicate programmes here; final database adjudication is unfinished.

Files retain site means, all component weights, 11-by-11 sampling covariance,
3-by-3 programme covariance, estimates and input hashes. Source details:
[Cieraad](../cieraad_sources/README.md), [Ceppi](../beech_sources/README.md),
[Guo raw source](https://doi.org/10.6084/m9.figshare.27273918.v1).
[Independent calculation](../search/leaf_damage_meta_review.md).
The main 108-effect catalogue is unchanged; this is an exploratory supplement.

## Direct support for the macro cohort

`plant_scope.csv` checks accepted-binomial exact matches against the fixed
58-candidate/305-comparator labels (361 distinct accepted taxa). None of the
four resolved experimental names matches that set. Aporosa remains unresolved
and is not coded as absent. This does not rule out future synonym matching.
European EIV MOIST is numeric for Fagus sylvatica (5.5) and Quercus robur (5.7),
non-numeric `x` for Vaccinium myrtillus, and unmatched for Psychotria asiatica.
These values are species indicators, not measured experimental-site moisture.
Consequently, the pooled result cannot directly establish or reject a
herbivory mechanism for the wet-decline candidates.

Reproduce the support audit and forest plot with
`python database/alan_experiments/leaf_damage_synthesis/check_scope_and_plot.py`.
The plot is exported as PNG, SVG and PDF under
`manuscript/biotic_synthesis/figures/leaf_damage_meta.*`.
