#!/usr/bin/env python3
"""Global ORIGINAL paired day/night fruit/seed field reproductive output vs prior
night-light at the ACTUAL author study site and experiment year (1993-2022).

This measures within-experiment night-versus-day reproductive ratios, but
historic corrected PCNL is still OBSERVATIONAL across different field studies.
No local lighting manipulation, no moth-specific pollen-vector identification,
no artificial extra positive/negative study cases, and no island-only claim.

Each original study×plant×metric×site×year×season is one source response
for descriptive WLS and study-cluster inference; studies receive equal total
weight independent of source number of repeated metrics/sites.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
from scipy.stats import t as student_t

LIGHT="pcnl_prior_mean_log1p"
SCHEMA="ALAN_DvN_reproductive_real_priorPCNL_source_join_v1"

def require(ok,why):
    if not bool(ok):raise ValueError("ALAN_TRUE_DIEL_REPRO_ALAN: "+str(why))

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def load_source(joined, receipt):
    r=json.loads(Path(receipt).read_text())
    require(r["schema"]==SCHEMA and
            r["original_source_pairs"]==301 and
            r["annual_rasters_MD5_source_verified"]==30 and
            sha(joined)==r["joined_file_sha256"] and
            r["true_source_geography_year_prior_light_only"],
            "Authenticated original paired experiment-year PCNL source altered")
    d=pd.read_csv(joined,dtype={"study_key":str,"species":str})
    require(len(d)==301 and d.source_pair_id.is_unique and
            d.year_exp.between(1993,2022).all() and
            d.night_share.between(0,1).all(),"Original source pair changed")
    valid=d.pcnl_historical_complete.astype(str).str.lower().eq("true")
    valid&=d.pcnl_prior_actual_sampled_years.eq(d.year_exp-1992)
    valid&=d.pcnl_prior_expected_years.eq(d.year_exp-1992)
    valid&=np.isfinite(pd.to_numeric(d[LIGHT],errors="coerce"))
    x=d.loc[valid].copy()
    x["lat2"]=x.lat.to_numpy(float)**2
    x["lon2"]=x.lon.to_numpy(float)**2
    x["exp_year"]=x.year_exp.astype(int)
    x["elevation_km"]=x.elevation_m.astype(float)/1000
    x["country"]=x.country.astype(str).replace({"NA":"unknown","":"unknown"})
    x["coordinate_source"]=x.coordinate_source.astype(str)
    x["geo20"]=((x.lat+90)//20).astype(int).astype(str).str.cat(
               ((x.lon+180)//20).astype(int).astype(str),sep="_")
    x["study_equal_weight"]=1/x.groupby("study_key",observed=True).study_key.transform("size")
    if x.empty:return x,r
    require(x.study_key.nunique()>=50 and x.species.nunique()>=60,
            "Insufficient independent external source studies to fit macro association")
    return x,r


def estimate(d, label, formula, field="night_share",two_way=True):
    if (len(d)<80 or d.study_key.nunique()<35 or
        d[LIGHT].std()<.02):
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_OR_ALAN_RANGE",
                "n_original_pairs":int(len(d)),"original_studies":int(d.study_key.nunique())}
    z=d.copy()
    z["study_equal_weight"]=1/z.groupby("study_key",observed=True).study_key.transform("size")
    if not np.isfinite(pd.to_numeric(z[field],errors="coerce")).all():
        return {"status":"NONFINITE_ORIGINAL_OUTCOME"}
    f=smf.wls(formula,data=z,weights=z.study_equal_weight).fit(
        cov_type="cluster",cov_kwds={"groups":z.study_key},use_t=False)
    rank=int(np.linalg.matrix_rank(f.model.exog))
    if LIGHT not in f.params or rank!=len(f.params):
        return {"status":"NOT_ESTIMABLE_DESIGN_OR_GEO_FIXED_EFFECT_RANK",
                "n_params":int(len(f.params)),"rank":rank}
    b=float(f.params[LIGHT])
    ci=f.conf_int().loc[LIGHT]
    x={
      "status":"EXPLORATORY_FIELD_STUDY_NIGHT_VS_DAY_REPRODUCTIVE_OUTPUT_PRIOR_PCNL",
      "model":label,"original_study_equal_weight":True,
      "response_original_within_experiment":field,
      "pcnl_prior_light_coefficient_fraction_per_mean_log1p":b,
      "95CI_study_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
      "p_nominal_study_cluster":float(f.pvalues[LIGHT]),
      "n_paired_original_experiment_locations_year_metrics":int(len(z)),
      "n_independent_original_studies":int(z.study_key.nunique()),
      "n_source_plant_species":int(z.species.nunique()),
      "n_distinct_original_sites":int(z[["lat","lon"]].drop_duplicates().shape[0]),
      "n_fixed_effect_parameters":int(len(f.params)),
      "geographic_control":formula,
      "night_only_response_is_not_specific_to_moths":True,
      "prior_light_at_site_is_not_experimental_light_manipulation":True,
    }
    if two_way:
        try:
            studies=pd.factorize(z.study_key)[0]
            species=pd.factorize(z.species)[0]
            v,_,_=cov_cluster_2groups(f,studies,species)
            beta_idx=list(f.params.index).index(LIGHT)
            variance=float(v[beta_idx,beta_idx])
            if np.isfinite(variance) and variance>0:
                se=float(np.sqrt(variance))
                df=min(z.study_key.nunique(),z.species.nunique())-1
                q=student_t.ppf(.975,df)
                x["two_way_independent_study_plant_cluster"]={
                    "SE":se,"CI95_student_t":[b-q*se,b+q*se],
                    "p_two_sided_t":float(2*student_t.sf(abs(b/se),df)),
                    "df":int(df)}
            else:x["two_way_independent_study_plant_cluster"]={"status":"NOT_ESTIMABLE_VARIANCE"}
        except Exception as e:
            x["two_way_independent_study_plant_cluster"]={"status":"NOT_ESTIMABLE","reason":str(e)[:120]}
    return x


def run(joined,receipt,output):
    output=Path(output)
    if output.exists():raise FileExistsError("Refuse replacing original audited result")
    d,src=load_source(joined,receipt)
    if d.empty:raise ValueError("No actual original point/corrected prior light overlap")
    methods={
      "original_unadjusted_across_study_sites":
        "night_share ~ "+LIGHT,
      "original_geography_experiment_year_metric_elevation":
        "night_share ~ "+LIGHT+
        " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)",
      "preexisting_geography_add_20degree_fixed_effects":
        "night_share ~ "+LIGHT+
        " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric) + C(geo20)"
    }
    fits={}
    for key,formula in methods.items():
        fits[key]=estimate(d,key,formula)
    # Actual DATE/LOCATION quality and author source native controls are not
    # picked based on sign or p-value; each is a separately labeled sensitivity.
    subgroups={
      "native_source_only":d.loc[d.plant_native_source.eq("native")],
      "author_original_coordinates_provided_in_text_only":
           d.loc[d.coordinate_source.eq("provided in text")],
      "non_multiyear_start_year_experiment_only":
           d.loc[d.year_end.isna()|d.year_end.eq(d.year_exp)],
      "fruit_set_only":d.loc[d.metric.eq("fruit set")],
      "seed_set_only":d.loc[d.metric.eq("seed set")],
      "hours_normalized_day_night_response":
           d.loc[d.night_share_per_hour.notna()],
    }
    for k,v in subgroups.items():
        if k=="hours_normalized_day_night_response":
            form=methods["original_geography_experiment_year_metric_elevation"].replace(
               "night_share ~ ","night_share_per_hour ~ ")
            fits[k]=estimate(v,k,form,field="night_share_per_hour")
        else:
            fits[k]=estimate(v,k,methods["original_geography_experiment_year_metric_elevation"])
    # Latitude, date and author-supplied elevation are PRE-LIGHT measured
    # independent confounders; absence of common on/off variation is a gate.
    rich=d.groupby("study_key",observed=True).agg(
          source_year_span=("year_exp",lambda z:int(z.max()-z.min())),
          site_count=("lat",lambda z:int(len(set(zip(z,d.loc[z.index,"lon"]))))),
          light_span=(LIGHT,lambda z:float(z.max()-z.min())),
          n_source_pairs=("source_pair_id","size"))
    support={
       "source_studies_with_2plus_field_experiment_site_years":int(rich.site_count.ge(2).sum()),
       "source_studies_with_2plus_experiment_years":int(rich.source_year_span.gt(0).sum()),
       "source_studies_with_2plus_light_values":int(rich.light_span.gt(.03).sum()),
       "source_studies_with_2plus_sites_and_light_variation":
           int((rich.site_count.ge(2)&rich.light_span.gt(.03)).sum()),
       "n_full_prior_original_pairs":len(d),
       "n_original_published_studies":int(d.study_key.nunique()),
       "n_original_field_sites":int(d[["lat","lon"]].drop_duplicates().shape[0]),
       "n_original_source_species":int(d.species.nunique()),
       "original_source_years":[int(d.year_exp.min()),int(d.year_exp.max())],
       "source_dark_PCNL_prior_mean_zero_count":int(d[LIGHT].eq(0).sum()),
       "source_light_prior_sd":float(d[LIGHT].std()),
       "source_night_share_median":float(d.night_share.median()),
       "source_night_share_mean_study_equal":float(np.average(d.night_share,weights=d.study_equal_weight)),
       "original_field_site_elevation_measured":int(d.elevation_m.notna().sum()),
       "source_author_night_hours_measured":int(d.night_share_per_hour.notna().sum()),
       "source_year_end_missing":int(d.experiment_year_end_unreported.sum()),
       "source_year_multi":int(d.year_end.gt(d.year_exp).sum()),
       "source_coordinate_quality_counts":{str(k):int(v) for k,v in d.coordinate_source.value_counts().items()},
       "source_native_pairs":int(d.plant_native_source.eq("native").sum()),
    }
    report={
      "schema":"ALAN_authentic_DvN_effective_day_vs_night_reproduction_301_field_siteyear_original_PRIOR_PCNL_v1",
      "status":"EXPLORATORY_ACTUAL_WITHIN_EXPERIMENT_DAYNIGHT_REPRO_AGAINST_HISTORICAL_LIGHT",
      "source":{
          "article_DOI":"10.1111/ele.70036",
          "author_original_treatment_blob":"39a82d966afcdde670603add1c881ace8791ffe9",
          "original_true_preexperiment_PCNL_years":[1992,2021],
          "original_source_pairs_before_light":301,
          "original_source_joined_SHA256":src["joined_file_sha256"],
          "original_30_annual_rasters_MD5_verified":True,
          "field_site_and_study_start_year_not_publication":True
      },
      "source_support":support,
      "exploratory_models":fits,
      "inferential_limits":{
          "day_night_exclusion_does_NOT_identity_moth_pollinators":True,
          "author_reproductive_outcome_is_measured_within_experiment":True,
          "illumination_not_experimentally_assigned_across_field_studies":True,
          "published_studies_nonrandomly_geographically_distributed":True,
          "global_source_not_only_island":True,
          "source_author_study_year_end_mostly_not_reported":True,
          "original_site_geo_precision_not_uniform":True,
          "fruit_and_seed_set_are_separate_metrics_within_paired_units":True,
          "study_and_species_source_nonindependence_clustered":True,
          "source_selected_by_geography_time_not_ALAN_or_night_share":True,
          "only_original_aggregate_result_exported_not_original_study_coordinates":True,
          "cannot_infer_historical_pollination_loss":True,
          "posthoc_meta_hypothesis_not_preregistered":True,
      }
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ACTUAL_TRUE_DAYNIGHT_FIELD_REPRO_PRIORPCNL_MACRO",
          json.dumps({"source":support,"models":{k:{
             "status":v.get("status"),
             "beta":v.get("pcnl_prior_light_coefficient_fraction_per_mean_log1p"),
             "ci":v.get("95CI_study_cluster"),"p":v.get("p_nominal_study_cluster"),
             "n":v.get("n_paired_original_experiment_locations_year_metrics"),
             "pub":v.get("n_independent_original_studies")}
            for k,v in fits.items()}},ensure_ascii=False),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("joined","receipt","output"):
        p.add_argument("--"+a,required=True,type=Path)
    args=p.parse_args()
    run(args.joined,args.receipt,args.output)
