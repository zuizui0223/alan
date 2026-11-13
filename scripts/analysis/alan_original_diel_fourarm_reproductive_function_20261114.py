#!/usr/bin/env python3
"""True four-arm day/night/24h-open/24h-excluded reproductive functional audit.

DATA ARE ORIGINAL, not a model of species-level moth labels:
- 301 SHA-authenticated field trial study×plant×site×year×season×fruit/seed
  pairs of true day-only vs night-only pollinator exclusion treatments.
- 1,652 original author treatment rows, Git blob verified. Extract open
  pollination AND complete-exclusion reproductive responses matching the
  SAME study, plant, site, field year/months, metric, and unit.
- Original 30 annual 1992–2021 MD5 PCNL sampled strictly BEFORE original
  field year; used only in a SECONDARY observational analysis.

PRIMARY: source 2x2 four-arm NONADDITIVITY of reproductive OUTPUT
  interaction = (Open + Closed - DayOnly - NightOnly) / (Open + Closed + DayOnly + NightOnly).
Negative = subadditivity of observed outputs on their ORIGINAL measurement scale.
This is not proof of biological redundancy: experimental bagging, wind/selfing,
saturation and nonrandom source-study geography also influence arm differences.

SECONDARY: opening both periods vs best of two single periods:
  open_increment=(Open - max(DayOnly,NightOnly)) /
                 (Open + max(DayOnly,NightOnly)).
Positive demonstrates numerical additional reproductive output of full-access
flowers in the source experiment, not verified pollen deposition or the
presence of a substitute moth/bee at the source field site.

Only aggregate results and diagnostics exported, no source study/plant/site
identities, field coordinates, copyright-restricted trait rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import binomtest

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, require, LIGHT,
)
from scripts.analysis.alan_flower_multitraits_reproductive_diel_transfer_20261113 import (
    original_traits,
)

AUTHOR_STUDY_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"
KEY=["study_ID","plant_species","treatment_effectiveness_metric",
     "coordinates_lat","coordinates_lon","year_start","year_end",
     "month_start","month_end","effectiveness_units"]
ARM=("day pollination","night pollination","open pollination","complete exclusion","hand pollination")
CORE_ARMS=ARM[:4]
SOURCE_N=(301,118,123)


def gitblob(path):
    data=Path(path).read_bytes()
    return hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()


def canonical_int(v):
    if pd.isna(v) or str(v) in ("","NA","nan","None"):return None
    f=float(v)
    require(np.isfinite(f) and f.is_integer(),"Original date not exact integer")
    return int(f)


def source_literal(v):
    # pandas default CSV reader converts the author's LITERAL "NA" units
    # and season labels to NaN. Restore those original source tokens
    # for exact byte-equivalent treatment group keys, never impute value.
    return "NA" if pd.isna(v) or str(v) in ("nan","None") else str(v).strip()


def canonical_key_from_raw(values):
    z=values.iloc[0]
    return (str(z.study_ID).strip()," ".join(str(z.plant_species).split()),
            str(z.treatment_effectiveness_metric),
            source_literal(z.effectiveness_units),
            round(float(z.coordinates_lat),7),
            round(float(z.coordinates_lon),7),
            canonical_int(z.year_start),canonical_int(z.year_end),
            source_literal(z.month_start),source_literal(z.month_end))


def canonical_key_from_preflight(row):
    return (str(row.study_key).strip()," ".join(str(row.species).split()),
            str(row.metric),source_literal(row.units),
            round(float(row.lat),7),round(float(row.lon),7),
            canonical_int(row.year_exp),canonical_int(row.year_end),
            source_literal(row.month_start),source_literal(row.month_end))


def original_mean(rows):
    if rows.empty:return None
    amount=pd.to_numeric(rows.effectiveness_value,errors="coerce")
    if amount.isna().any() or (~np.isfinite(amount)).any() or amount.lt(0).any():
        return None
    n=pd.to_numeric(rows.sample_size,errors="coerce")
    if n.notna().all() and np.isfinite(n).all() and n.gt(0).all():
        return float(np.average(amount.to_numpy(float),weights=n.to_numpy(float)))
    return float(amount.mean())


def exact_source_arms(original:pd.DataFrame,source:pd.DataFrame):
    """Fail closed if ANY of 301 original paired D+N source experiments
    cannot be reconstructed from literal original author treatment rows."""
    require(gitblob(original)==AUTHOR_STUDY_BLOB,
            "Original Kendall Nicholson public study table Git blob altered")
    r=pd.read_csv(original,dtype=str,keep_default_na=False)
    require(len(r)==1652 and r.study_ID.nunique()==140,
            "Published 1652×140 study treatment source changed")
    require(set(KEY+["treatment_condition","sample_size","effectiveness_value","exclude"])<=set(r),
            "Original published intervention treatment field absent")
    require(r.exclude.eq("n").all(),"Originally excluded study treatments must not be assumed included")
    bins={}
    for _,a in r.loc[r.treatment_effectiveness_metric.isin(["fruit set","seed set"])].groupby(
        KEY,dropna=False,observed=True,sort=False):
        first=a.iloc[0]
        try:key=canonical_key_from_raw(a)
        except (ValueError,TypeError):continue
        if not (1993<=key[6]<=2022 and
                np.isfinite(key[4]) and np.isfinite(key[5]) and
                abs(key[4])<=90 and abs(key[5])<=180):
            continue
        require(key not in bins,"Canonical source identifier merged two distinct original groups")
        bins[key]=a.copy()
    hit=set()
    outcome=[]
    missing=[]
    bad_daynight=0
    present_by_arm={t:0 for t in ARM}
    for row in source.itertuples(index=False):
        k=canonical_key_from_preflight(row)
        group=bins.get(k)
        if group is None:
            missing.append(True)
            continue
        require(k not in hit,"Duplicate source study×plant×site/year/season pair")
        hit.add(k)
        average={name:original_mean(group[group.treatment_condition.eq(name)])
                 for name in ARM}
        require(average["day pollination"] is not None and
                average["night pollination"] is not None,
                "Original site/time day or night source claim incomplete")
        # Night and day means already SHA-validated independently in parent,
        # so all 301 originals must match exactly and no "new" pair is invented.
        if not (np.isclose(average["day pollination"],float(row.day_reproduction),
                           atol=1e-9,rtol=1e-10) and
                np.isclose(average["night pollination"],float(row.night_reproduction),
                           atol=1e-9,rtol=1e-10)):
            bad_daynight+=1
        for name in ARM:
            if (group.treatment_condition.eq(name)).any():
                present_by_arm[name]+=1
        out={
            "source_pair_id":int(row.source_pair_id),
            "study_key":str(row.study_key),
            "species":str(row.species),
            "metric":str(row.metric),
            "year_exp":int(row.year_exp),
            "lat":float(row.lat),"lon":float(row.lon),
            "open":average["open pollination"],
            "closed":average["complete exclusion"],
            "hand":average["hand pollination"],
            "day":average["day pollination"],
            "night":average["night pollination"],
        }
        outcome.append(out)
    print("AUTHENTIC_FIELD_FOURARM_PAIR_KEY_AUDIT",json.dumps({
        "original_301_expected":SOURCE_N[0],"matched_source_pairs":len(outcome),
        "source_pair_groups_not_identically_reconstructed":len(missing),
        "source_day_night_mean_differences":bad_daynight,
        "original_valid_group_keys":len(bins)}),flush=True)
    require(len(missing)==0 and len(outcome)==SOURCE_N[0] and bad_daynight==0,
            "Original signed day/night source cannot be reconstructed from author interventions")
    for original_col,counted in (("open_pollination_arms_present","open pollination"),
                                 ("hand_pollination_arms_present","hand pollination")):
        require(int(source[original_col].sum())==present_by_arm[counted],
                "Original preflight intervention flag mismatches study source arms")
    require(present_by_arm["open pollination"]==264 and
            present_by_arm["hand pollination"]==126,
            "Original study site-year open/hand outcome coverage changed")
    all_arms=pd.DataFrame(outcome)
    q=source.merge(all_arms.drop(columns=["study_key","species","metric","year_exp","lat","lon"]),
                   on="source_pair_id",validate="one_to_one")
    require(len(q)==SOURCE_N[0],"Original real field source intervention expansion detected")
    documented_all=q[["open","closed","day","night"]].notna().all(axis=1)
    three=q[["open","day","night"]].notna().all(axis=1)
    hand=q[["open","hand"]].notna().all(axis=1)
    return q,{
        "original_source_treatment_rows":len(r),
        "original_source_studies":int(r.study_ID.nunique()),
        "original_reproductive_daynight_pairs":SOURCE_N[0],
        "original_open_24h_access_pairs_documented":present_by_arm["open pollination"],
        "original_complete_exclusion_24h_pairs_documented":present_by_arm["complete exclusion"],
        "original_hand_pollination_pairs_documented":present_by_arm["hand pollination"],
        "daynight_open_full_three_arms_matched":int(three.sum()),
        "daynight_open_closed_FULL_2x2_matched":int(documented_all.sum()),
        "daynight_open_hand_original_pairs":int(hand.sum()),
        "study_count_with_open_both":int(q.loc[three,"study_key"].nunique()),
        "study_count_full_2x2":int(q.loc[documented_all,"study_key"].nunique()),
        "species_count_full_2x2":int(q.loc[documented_all,"species"].nunique()),
        "study_daynight_original_mean_disagreements":bad_daynight,
        "original_author_study_gitblob":AUTHOR_STUDY_BLOB,
        "open_and_closed_are_ORIGINAL_author_defined_treatment_categories":True,
    }


def qualify(q:pd.DataFrame):
    """A genuine four-arm source with matched reproductive metric and >0
    observed original total outcome. Never divide by zero or invent a
    fifth intervention in an unobserved source study."""
    a=q[q[["open","closed","day","night"]].notna().all(axis=1)].copy()
    for col in ("open","closed","day","night"):
        require(np.isfinite(a[col].to_numpy(float)).all() and (a[col]>=0).all(),
                "Non-finite or negative author fruit/seed reproductive response")
    denom=a[["open","closed","day","night"]].sum(axis=1)
    zero=int(denom.eq(0).sum())
    a=a.loc[denom.gt(0)].copy()
    n=a.open.to_numpy(float)
    d=a.day.to_numpy(float)
    z=a.night.to_numpy(float)
    c=a.closed.to_numpy(float)
    a["open_relative_to_best"]=np.divide(n-np.maximum(d,z),n+np.maximum(d,z),
                                          out=np.full(len(a),np.nan),where=(n+np.maximum(d,z))>0)
    a["two_window_nonadditivity"]=(n+c-d-z)/(n+c+d+z)
    a["closed_reproduction_fraction"]=c/(n+c+d+z)
    a["open_vs_closed"]=np.divide(n-c,n+c,out=np.full(len(a),np.nan),where=(n+c)>0)
    a["two_window_source_scale_total"]=n+c+d+z
    require(a["two_window_nonadditivity"].between(-1.,1.).all(),"Original bounded nonadditive index invalid")
    require(a.open_relative_to_best.dropna().between(-1,1).all(),
            "Original 24h vs best single access outcome invalid")
    # If only the closed outcome is positive, the open/max-single ratio
    # is 0/0 and must remain unknown, not an invented zero effect.
    require((a.study_key.nunique()>=50 and a.species.nunique()>=60),
            "Too few independent original fully measured four-arm studies/taxa")
    return a,{
        "all_original_four_arm_available":int(len(a)+zero),
        "total_outcome_all_four_arms_zero_EXCLUDED":zero,
        "valid_four_arm_pairs":len(a),
        "independent_studies":int(a.study_key.nunique()),
        "independent_plant_species":int(a.species.nunique()),
        "original_outcome_metric_counts":{str(k):int(v) for k,v in a.metric.value_counts().items()},
        "source_site_year_native_bloom_not_confirmed":True,
    }


def descriptive_study_sign(a:pd.DataFrame,metric:str):
    """Sign inference at original publication grain, not duplicated field
    seasons from a single paper treated as independent."""
    require(metric in ("open_relative_to_best","two_window_nonadditivity","open_vs_closed"),
            "Unsupported source intervention contrast")
    vals=a[["study_key","species",metric]].dropna()
    # Species is not forced unique per publication; report medians by
    # study, and independently by species for dual source leverage.
    study=vals.groupby("study_key")[metric].median()
    species=vals.groupby("species")[metric].median()
    n_neg=int(study.lt(-1e-12).sum());n_pos=int(study.gt(1e-12).sum())
    n_tie=len(study)-n_neg-n_pos
    if len(study)<30 or len(species)<35:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_BIOLOGICAL_STUDY_SUPPORT",
                "original_studies":int(len(study)),"species":int(len(species))}
    test=binomtest(n_neg,n_neg+n_pos,.5,alternative="two-sided") if n_neg+n_pos else None
    # Publication-level descriptive bootstrap on original publications;
    # not phylogenetically adjusted, and only published study populations.
    rng=np.random.default_rng(20261114)
    v=study.to_numpy(float)
    indices=rng.integers(0,len(v),(3000,len(v)))
    med=np.median(v[indices],axis=1)
    return {
      "status":"ACTUAL_ORIGINAL_FOUR_ARM_STUDY_GRAIN_DESCRIPTIVE",
      "index":metric,
      "source_original_comparisons":int(len(vals)),
      "source_studies":len(study),
      "source_species":len(species),
      "median_all_original_reproductive_pairs":float(vals[metric].median()),
      "median_original_independent_study_effects":float(study.median()),
      "median_independent_species_effects":float(species.median()),
      "mean_original_independent_study_effects":float(study.mean()),
      "original_independent_study_median_negative":n_neg,
      "original_independent_study_median_positive":n_pos,
      "original_independent_study_median_exactly_zero":n_tie,
      "source_study_sign_test_two_sided_p_nominal":float(test.pvalue) if test else None,
      "study_level_bootstrap_median_95pct_DESCRIPTIVE_interval":[float(np.quantile(med,.025)),
                                                                   float(np.quantile(med,.975))],
      "study_bootstrap_original_publication_groups":3000,
      "not_3000_independent_replications":True,
    }


def original_observational_pcNL(a:pd.DataFrame,metric):
    """Secondary exposure association only; true light not experimentally
    assigned. Study-equal weights, independent publication clusters."""
    x=a.loc[a[metric].notna()].copy()
    if len(x)<110 or x.study_key.nunique()<55 or x[LIGHT].std()<.02:
        return {"status":"NOT_ESTIMABLE_TRUE_STUDY_LIGHT_SUPPORT",
                "observations":len(x),"studies":int(x.study_key.nunique())}
    x["exp_year"]=x.year_exp.astype(int)
    x["lat2"]=x.lat**2;x["lon2"]=x.lon**2
    x["elevation_km"]=x.elevation_m.astype(float)/1000
    x["study_equal_weight"]=1/x.groupby("study_key").study_key.transform("size")
    formula=(metric+" ~ "+LIGHT+" + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")
    # Metric name also 'metric' for reproductive fruit/seed treatment. Use
    # actual response alias to avoid collision with C(metric) categorical.
    formula=("source_response ~ "+LIGHT+
             " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")
    x["source_response"]=x[metric].to_numpy(float)
    fit=smf.wls(formula,data=x,weights=x.study_equal_weight).fit(
        cov_type="cluster",cov_kwds={"groups":x.study_key},use_t=False)
    if np.linalg.matrix_rank(fit.model.exog)!=len(fit.params):
        return {"status":"NOT_ESTIMABLE_FIXED_EFFECT_RANK",
                "studies":int(x.study_key.nunique())}
    ci=fit.conf_int().loc[LIGHT]
    return {"status":"EXPLORATORY_ORIGINAL_WITHIN_FIELD_REPRODUCTIVE_CONTRAST_ACROSS_STUDY_LIGHT",
            "original_response":metric,
            "original_study_pairs":len(x),
            "published_studies":int(x.study_key.nunique()),
            "plant_species":int(x.species.nunique()),
            "study_equal_weight":True,
            "beta_fraction_per_mean_log1p_prior_PCNL":float(fit.params[LIGHT]),
            "95CI_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
            "p_nominal":float(fit.pvalues[LIGHT]),
            "original_year_geography_elevation_metric_adjusted":True,
            "not_randomized_at_flower_illumination":True,
            "no_nocturnal_moth_specific_effect":True}


def full_analysis(joined,receipt,author,plant_traits,out):
    out=Path(out)
    if out.exists():raise FileExistsError("Do not overwrite original published-source audit")
    source,original_receipt=load_source(joined,receipt)
    require((len(source),source.study_key.nunique(),source.species.nunique())==SOURCE_N,
            "Original 301 true site/year day/night experiments changed")
    allsource,coverage=exact_source_arms(author,source)
    a,eligible=qualify(allsource)
    trait=original_traits(plant_traits)
    b=a.merge(trait[["species","bloom_period_simple","plant_family"]],
              on="species",validate="many_to_one",how="left")
    group_support={
        "documented_day_bloom_fourarm_pairs":int(b.bloom_period_simple.eq("day").sum()),
        "documented_night_bloom_fourarm_pairs":int(b.bloom_period_simple.eq("night").sum()),
        "native_field_source_fourarm_pairs":int(b.plant_native_source.eq("native").sum()),
        "original_reported_day_night_hours_fourarm_pairs":int(b.night_share_per_hour.notna().sum()),
        "original_dark_satellite_prior_PCNL_pairs":int(b[LIGHT].eq(0).sum()),
        "original_nonzero_satellite_prior_PCNL_pairs":int(b[LIGHT].gt(0).sum()),
        "metric_source_counts":eligible["original_outcome_metric_counts"],
    }
    primary={
        name:descriptive_study_sign(a,name)
        for name in ("two_window_nonadditivity","open_relative_to_best","open_vs_closed")
    }
    stratified={}
    for name,subset in {
       "source_native_plant_only":a.loc[a.plant_native_source.eq("native")],
       "reported_exposure_hours_only":a.loc[a.night_share_per_hour.notna()],
       "fruit_set_only":a.loc[a.metric.eq("fruit set")],
       "seed_set_only":a.loc[a.metric.eq("seed set")],
    }.items():
        stratified[name]={
           "original_comparisons":len(subset),
           "source_studies":int(subset.study_key.nunique()),
           "source_species":int(subset.species.nunique()),
           "nonadditivity":descriptive_study_sign(subset,"two_window_nonadditivity"),
           "open_vs_best":descriptive_study_sign(subset,"open_relative_to_best")
        }
    output={
       "schema":"ALAN_real301_diel_fourarm_field_experiment_reproductive_insurance_PCNL_source_v1",
       "status":"SOURCE_AUTHENTIC_ORIGINAL_DAY_NIGHT_OPEN_CLOSED_FRUIT_SEED_CONTRAST",
       "source":{
          "author_DOI":"10.1111/ele.70036",
          "author_original_study_Gitblob":AUTHOR_STUDY_BLOB,
          "true_daynight_open_closed_site_year_source_PCNLP_SHA256":original_receipt["joined_file_sha256"],
          "original_full_301_pairs":SOURCE_N[0],
          "original_independent_studies":SOURCE_N[1],
          "original_plant_species":SOURCE_N[2],
          "previous_light_annual_rasters_MD5_verified":[1992,2021],
       },
       "arm_source_coverage":coverage,
       "four_arm_positive_total_source_support":eligible,
       "fourarm_other_study_metadata":group_support,
       "primary_true_experimental_four_arm_reproductive_contrasts":primary,
       "original_functional_sensitivity_by_protocol":stratified,
       "secondary_prior_satellite_light_associations":{
           name:original_observational_pcNL(a,name) for name in
           ("two_window_nonadditivity","open_relative_to_best")
       },
       "biological_interpretation_guardrails":{
          "open_can_include_unmeasured_wind_or_visit_behaviour":True,
          "closed_can_suppress_wind_and_alter_microclimate":True,
          "fruit_seed_maximum_causes_saturation_of_additivity":True,
          "original_day_night_exclusion_does_not_ID_moth_taxa":True,
          "satellite_prior_light_not_manipulated_illumination":True,
          "publication_source_selection_nonrandom":True,
          "study_median_sign_tests_are_not_species_phylogeny_models":True,
          "separate_experimental_open_vs_single_comparisons_already_in_parent_meta_analysis":True,
          "no_true_local_pollinator_shift_or_autonomous_selfing_compensation_proven":True,
          "no_restricted_original_plant_study_site_data_exported":True,
       },
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_GENUINE_FOURARM_DIEL_FUNCTION",
          json.dumps({"coverage":coverage,"eligible":eligible,
                      "main":primary,"strata":stratified,
                      "pcnl":output["secondary_prior_satellite_light_associations"]},
                      ensure_ascii=False),flush=True)
    return output


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    for n in ("joined","receipt","author","trait","out"):
        a.add_argument("--"+n,type=Path,required=True)
    z=a.parse_args()
    full_analysis(z.joined,z.receipt,z.author,z.trait,z.out)
