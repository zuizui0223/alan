#!/usr/bin/env python3
"""Exploratory measured pollen-limitation vs actual pre-experiment ALAN.

This is a publication-clustered study-equal-weight observational comparison.
The exposure is real MD5-audited 1992–2014 PCNL sampled at GloPL experiment
coordinates in years STRICTLY before year_of_study (not publication year).
Pollen limitation is the published hand-pollen vs natural reproduction log
response ratio. It is NOT moth-only pollination, NPMS plant cover, a measured
ground lux exposure, or evidence of causal ALAN-driven pollen deficiency.

The model and support thresholds were set before inspecting the joined PL
response. Do not tune covariates or spatial blocks to find significance.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

EXPOSURE="pcnl_prior_mean_log1p"
MODEL_SCHEMA="ALAN_real_original_GloPL_pollensupplementation_vs_priorPCNL_exploratory_v1"
SOURCE_SCHEMA="ALAN_real_GloPL_PCNLPrior1992_2014_source_join_v1"
METH_COLS=["PL_Effect_Size_Type1","PL_Effect_Size_Type2",
           "Constant_added","Level_of_Supplementation"]


def require(ok,message):
    if not ok:raise ValueError("ALAN_GLOPL_PL_ASSOCIATION: "+str(message))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()


def validate(joined,receipt):
    report=json.loads(Path(receipt).read_text())
    require(report["schema"]==SOURCE_SCHEMA and
            report["source_experiment_rows"]==1419 and
            report["original_PCNLP_years_verified"]==23 and
            report["GloPL_PL_outcome_not_used_to_choose_source_cohort"] and
            report["historical_light_never_future_experiment"] and
            sha(joined)==report["output_sha256"],
            "Original source light+PL cohort must be verified")
    x=pd.read_csv(joined,dtype={"study_key":str,"species":str},
                  low_memory=False)
    required={"original_row_id","study_key","species","year_exp",
              "lat","lon","PL_Effect_Size",EXPOSURE,
              "pcnl_historical_complete","pcnl_pre_experiment_expected_years",
              "pcnl_pre_experiment_observed_years",*METH_COLS}
    require(len(x)==1419 and required.issubset(x.columns) and
            x.original_row_id.is_unique,"GloPL original experimental data altered")
    require((x.pcnl_pre_experiment_expected_years==x.year_exp-1992).all() and
            (x.pcnl_pre_experiment_observed_years <=
             x.pcnl_pre_experiment_expected_years).all(),
            "Future or double-counted light observations")
    return x,report


def gate_source(x):
    valid=x.pcnl_historical_complete.astype(str).str.lower().eq("true")
    valid&=x.pcnl_pre_experiment_observed_years.eq(
        x.pcnl_pre_experiment_expected_years)
    valid&=np.isfinite(pd.to_numeric(x[EXPOSURE],errors="coerce"))
    # Source cohort is independent of PL effect size: outcome completeness is
    # explicitly counted only after source-year and exposure support.
    return x.loc[valid].copy()


def fit(g, formula, study_clusters):
    if len(g)<200 or g.study_key.nunique()<50 or g[EXPOSURE].std()<.03:
        return {"status":"insufficient_source_support",
                "n_rows":len(g),
                "n_studies":g.study_key.nunique(),
                "light_sd":float(g[EXPOSURE].std())}
    m=smf.wls(formula,data=g,weights=g.publication_equal_weight).fit(
          cov_type="cluster", cov_kwds={"groups":g.study_key},use_t=False)
    name=EXPOSURE
    require(name in m.params and np.isfinite(m.bse[name]) and
            len(m.params) < len(g)/2,
            "Effect not identifiable from model or parameters overfit")
    b=float(m.params[name]);se=float(m.bse[name])
    ci=m.conf_int().loc[name]
    return {
      "status":"estimated_exploratory_not_causal",
      "beta_GloPL_log_pollen_supplementation_ratio_per_meanLogPCNL":b,
      "publication_cluster_robust_SE":se,
      "95CI":[float(ci.iloc[0]),float(ci.iloc[1])],
      "nominal_two_sided_p":float(m.pvalues[name]),
      "n_rows":int(m.nobs),
      "n_unique_publications":int(g.study_key.nunique()),
      "n_unique_coordinates":int(g[["lat","lon"]].drop_duplicates().shape[0]),
      "n_unique_taxon_names":int(g.species.nunique()),
      "n_predictors":len(m.params),
      "model_R2_descriptive":float(m.rsquared),
      "formula":formula,
      "response_label":"published pollen-supplementation experiment log response ratio",
      "not_effect_of_ALAN_on_NPMS_plant_occupancy_or_nocturnal_moths":True
    }


def run(joined,receipt,out):
    x,source=validate(joined,receipt)
    g=gate_source(x)
    y=pd.to_numeric(g.PL_Effect_Size,errors="coerce")
    finite=np.isfinite(y)
    n_nonfinite=int((~finite).sum())
    g=g.loc[finite].copy()
    g["PL_log_ratio"]=y.loc[g.index]
    g[EXPOSURE]=pd.to_numeric(g[EXPOSURE])
    g["exp_year"]=g.year_exp.astype(int)
    g["lat2"]=g.lat.astype(float)**2
    g["lon2"]=g.lon.astype(float)**2
    g["lat"]=g.lat.astype(float)
    g["lon"]=g.lon.astype(float)
    for col in METH_COLS:
        # Actual categorical level preserved, not assigned arbitrary biology.
        g[col]=g[col].astype(str).fillna("unknown").replace({"nan":"unknown","":"unknown"})
    total=g.groupby("study_key",observed=True).size()
    g["publication_equal_weight"]=1.0/g.study_key.map(total).to_numpy(float)
    g["geo10"]=((g.lat+90)//10).astype(int).astype(str).str.cat(
               ((g.lon+180)//10).astype(int).astype(str),sep="_")
    light_q=g[EXPOSURE].quantile([.1,.25,.5,.75,.9]).to_dict() if len(g) else {}
    source_support={
      "rows_original":len(x),
      "rows_full_priorPCNL":int(gate_source(x).shape[0]),
      "rows_finite_experiment_PL_and_full_PCNL":len(g),
      "rows_missing_PL_with_complete_PCNL":n_nonfinite,
      "n_unique_pub":int(g.study_key.nunique()),
      "n_unique_original_GloPL_sites":int(g[["lat","lon"]].drop_duplicates().shape[0]),
      "n_10deg_spatial_blocks":int(g.geo10.nunique()),
      "n_unique_years":int(g.exp_year.nunique()),
      "year_min":int(g.exp_year.min()) if len(g) else None,
      "year_max":int(g.exp_year.max()) if len(g) else None,
      "fraction_zero_historical_PCNL_mean":float(g[EXPOSURE].eq(0).mean()) if len(g) else None,
      "mean_prior_light_log1p":float(g[EXPOSURE].mean()) if len(g) else None,
      "sd_prior_light_log1p":float(g[EXPOSURE].std()) if len(g) else None,
      "quartiles_prior_light":{str(k):float(v) for k,v in light_q.items()},
      "PL_effect_mean_descriptive":float(g.PL_log_ratio.mean()) if len(g) else None,
      "PL_effect_sd_descriptive":float(g.PL_log_ratio.std()) if len(g) else None,
    }
    methods=" + ".join(f"C({c})" for c in METH_COLS)
    base="PL_log_ratio ~ "+EXPOSURE+" + exp_year + lat + lat2 + lon + lon2 + "+methods
    specs={
      "global_study_equal_method_and_geography_adjusted":base,
      "global_add_10deg_geographical_block_FE":base+" + C(geo10)"
    }
    result={}
    if len(g) and (g[EXPOSURE].quantile(.75)-g[EXPOSURE].quantile(.25) > .02):
        for label,formula in specs.items():
            try:
                result[label]=fit(g,formula,"study_key")
            except (np.linalg.LinAlgError,ValueError) as error:
                result[label]={"status":"non_estimable",
                               "reason":str(error)[:250]}
    else:
        result={"global_study_equal_method_and_geography_adjusted":{
            "status":"insufficient_actual_historical_light_variation"}}
    report={
      "schema":MODEL_SCHEMA,
      "original_outcome_blind_source_sha256":source["original_outcome_blind_source_sha256"],
      "original_joined_actual_light_source_sha256":source["output_sha256"],
      "source_support":source_support,
      "model_results":result,
      "analysis_is_observational_and_exploratory":True,
      "outcome_is_supplementation_advantage_not_actual_nocturnal_moth_dependence":True,
      "NPMS_plant_cover_not_directly_linked_to_GloPL_experiment":True,
      "GloPL_publication_cluster_NOT_plant_populations":True,
      "historic_light_experiment_year_not_PUBLICATION_YEAR":True,
      "spatial_temporal_landuse_and_study_selection_confounding_unresolved":True,
      "not_a_causal_mediation_analysis":True,
    }
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of independently rechecked source report")
    out.mkdir(parents=True)
    (out/"real_GloPL_experimental_PL_vs_priorPCNL_exploratory_result.json").write_text(
         json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_GLOPL_PL_VS_HISTORIC_ALAN_EXPLORATORY",
          json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("joined","receipt","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.joined,a.receipt,a.out)


if __name__=="__main__":main()
