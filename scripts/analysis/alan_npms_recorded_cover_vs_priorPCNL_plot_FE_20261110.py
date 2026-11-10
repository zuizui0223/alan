#!/usr/bin/env python3
"""Real NPMS repeated-positive floral cover versus historical prior-year light.

Only original 2015–25 NPMS site×plant×real-survey-ID×season×year positive
cover records with the source CRC-verified DOMIN -> original [0,1] midpoint,
and actual MD5-audited 2014–24 PCNL year-minus-one source exposures.
Primary conditional-on-positive-recorded-plant cover (percent midpoint)
~ historical prior-year log1p(PCNL) + original day-of-year, with
plot×species×original-surveyID×season fixed effects AND species×calendar-year
fixed effects; cluster-robust uncertainty at REAL NPMS plot, equal-plot
weights. Nonrecord of species is UNKNOWN, not abundance=0 or extinction.
A species-specific independent GloPL PL effect was not measured here.
No causal nightlight/flower moth pollen/seed conclusions are made.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
from scipy.stats import t

SCHEMA="ALAN_original_NPMS_recorded_cover_vs_yearminus1_PCNL_withinplot_v1"
AUDIT_SCHEMA="ALAN_actual_NPMS_sameplot_survey_season_cover_priorPCNL_support_v1"
SOURCE_PCNl_SHA="5ecae6a9b62de636c8e68fd694c69698eec154e5a8a0ba641beb5a7e6061e427"
EARLIEST=2015;LATEST=2025;MIN_SOURCE_LIGHT_RANGE=.02
SOURCE_COLUMNS=["location_id","survey_method_id","season","species","year",
                "median_cover","median_source_day_of_year",
                "n_visits_with_cover","source_domin_encodings",
                "source_record_status","original_annual_priorPCNL",
                "source_prior_unique_values"]

def require(cond,msg):
    if not bool(cond):raise ValueError("ALAN_NPMS_COVER_TRUE: "+str(msg))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as inp:
        for ch in iter(lambda:inp.read(4<<20),b""):h.update(ch)
    return h.hexdigest()

def read_original(csv,receipt):
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==AUDIT_SCHEMA and
            r.get("original_16960_plotyear_light_SHA256")==SOURCE_PCNl_SHA and
            r.get("source_CRC_pinned_DOMIN_midPoint_unit")==
              "fraction of cover [0,1], not percent [0,100]" and
            r.get("source_records_have_actual_doy_for_within_season_phenology_adjustment") and
            r.get("source_original_positive_cover_table_SHA256")==sha(csv),
            "No real signed original NPMS DOMIN/fixed-plot/year-source audit")
    x=pd.read_csv(csv,dtype={"location_id":str,"survey_method_id":str,
                             "species":str,"season":str},low_memory=False)
    require(set(SOURCE_COLUMNS).issubset(x.columns) and
            not x.duplicated(["location_id","survey_method_id",
                              "season","species","year"]).any() and
            x.year.between(EARLIEST,LATEST).all() and
            x.season.isin(["spring","summer"]).all(),
            "An original NPMS source-positive plot×plant×year is duplicated/unknown")
    require(x.median_cover.between(0,1).all() and
            x.median_source_day_of_year.between(1,366).all() and
            x.original_annual_priorPCNL.ge(0).all() and
            (x.source_prior_unique_values==1).all() and
            x.n_visits_with_cover.ge(1).all(),
            "Original published 20 cover codes, DOY or pre-year PCNL corrupted")
    require(x.location_id.nunique()<=1203,
            "Cannot have more NPMS original light-gradient plots than frozen source")
    x["panel_key"]=(
        x.location_id+"|"+x.survey_method_id+"|"+x.season+"|"+x.species)
    x["species_year_key"]=x.species+"|"+x.year.astype(str)
    x["recorded_cover_pp"]=100.*x.median_cover.astype(float)
    x["recorded_cover_log"]=np.log1p(x.recorded_cover_pp)
    x["site_light"]=x.original_annual_priorPCNL.astype(float)
    x["doy_30"]=(x.median_source_day_of_year.astype(float)-180.)/30.
    x["log_recorded_visits"]=np.log1p(x.n_visits_with_cover.astype(float))
    return x,r

def choose(x,source_min_range=.02, same_domin_encoding=True,
           source_status_stable=False, max_date_deviation=None,
           min_num_years=2, only_survey_id=None):
    z=x.copy()
    if only_survey_id is not None:
        z=z.loc[z.survey_method_id.eq(str(only_survey_id))].copy()
    if max_date_deviation is not None:
        center=z.groupby("panel_key").median_source_day_of_year.transform("median")
        z=z.loc[(z.median_source_day_of_year-center).abs().le(max_date_deviation)].copy()
    if z.empty:return z,{"status":"no_source_positive_cover_rows"}
    groups=z.groupby("panel_key",sort=False).agg(
        source_true_observation_years=("year","nunique"),
        original_radiance_levels=("site_light","nunique"),
        min_light=("site_light","min"),
        max_light=("site_light","max"),
        original_domin_scheme_categories=("source_domin_encodings","nunique"),
        source_status_categories=("source_record_status","nunique"))
    valid=(groups.source_true_observation_years>=min_num_years)&(
        groups.original_radiance_levels>=2)&(
        groups.max_light-groups.min_light>=source_min_range)
    if same_domin_encoding:valid&=groups.original_domin_scheme_categories.eq(1)
    if source_status_stable:valid&=groups.source_status_categories.eq(1)
    eligible=z.loc[z.panel_key.isin(groups.index[valid])].copy()
    report={
       "status":"source_supported_for_same_plot_source_positive_cover",
       "source_positive_cover_rowyears":len(eligible),
       "source_true_independent_NPMS_plots":int(eligible.location_id.nunique()),
       "source_true_source_panel_plot_species_method_season_groups":int(valid.sum()),
       "distinct_positive_recorded_plant_species":int(eligible.species.nunique()),
       "distinct_source_years":sorted(eligible.year.unique().astype(int).tolist()),
       "original_observer_method_category_counts":{
           str(k):int(v) for k,v in eligible.survey_method_id.value_counts().items()},
       "source_light_range_cutoff_log1p":source_min_range,
       "require_original_domin_encoding_same_across_all_years":same_domin_encoding,
       "require_original_record_status_same_across_years":source_status_stable,
       "max_source_visit_DOY_difference_from_plot_plant_season_median":
          max_date_deviation,
       "min_two_positive_recording_years_not_plant_occurrence_absence":True}
    return eligible,report

def weighted_panel_partial_out(Z,groups,weight,tol=1e-8,maxit=2000):
    """Project original source y/X against *two sets* of true source FEs.

    Sparse LSMR on weighted photo-source group/year dummy variables, not a
    dense or pseudoreplicated original 17k-panel fixed-effect matrix.
    Zero-contributor original singleton groups partial out automatically.
    """
    n,k=Z.shape
    require(n>4 and len(groups)==2 and np.all(np.isfinite(Z)) and
            np.all(weight>0) and np.isfinite(weight).all(),
            "Original within-study design or weights invalid")
    codes=[pd.factorize(v,sort=True)[0].astype(np.int32) for v in groups]
    levels=[int(c.max()+1) for c in codes]
    rr=np.tile(np.arange(n,dtype=np.int32),len(codes))
    cc=np.concatenate([code+sum(levels[:i]) for i,code in enumerate(codes)])
    ff=csr_matrix((np.ones(len(rr),dtype=float),(rr,cc)),
                  shape=(n,sum(levels)),dtype=float)
    sqrt=np.sqrt(weight)
    fw=ff.multiply(sqrt[:,None]).tocsr()
    residual=np.zeros_like(Z,dtype=float)
    lsmr_evidence=[]
    for col in range(k):
        q=Z[:,col]*sqrt
        fit=lsmr(fw,q,atol=tol,btol=tol,maxiter=maxit)
        residual[:,col]=Z[:,col]-ff@fit[0]
        lsmr_evidence.append({
            "source_column_index":col,
            "source_sparse_LSMR_termination_code":int(fit[1]),
            "source_sparse_LSMR_iterations":int(fit[2]),
            "source_sparse_LSMR_relative_norm":float(fit[3]/max(np.linalg.norm(q),1e-12))})
        require(int(fit[1]) not in (0,3,6,7) and int(fit[2])<maxit,
                "Genuine within-plot/source year FE absorption not converged")
    require(np.isfinite(residual).all(),"Nonfinite source-positive fixed effects")
    return residual,{"original_FE_panel_groups":levels[0],
                     "original_FE_species_year_groups":levels[1],
                     "original_FE_dummy_nonzeros":int(ff.nnz),
                     "source_lsmr_numerical_diagnostics":lsmr_evidence}

def fit_source(d,y="recorded_cover_pp",equal="plot_equal",
               yearFE="species_year",include_doy=True):
    N=len(d)
    G=int(d.location_id.nunique())
    if N<400 or G<35 or d.year.nunique()<3:
        return {"status":"source_insufficient_repeat_positive_cover", 
            "original_positive_cover_rows":N,"source_true_plots":G}
    require(y in ("recorded_cover_pp","recorded_cover_log","log_recorded_visits"),
            "Unknown outcome would be post-hoc")
    if equal=="plot_equal":
        w=1/d.groupby("location_id").location_id.transform("size").to_numpy(float)
    elif equal=="plot_species_equal":
        w=1/d.groupby("panel_key").panel_key.transform("size").to_numpy(float)
    elif equal=="record_equal":
        w=np.ones(N)
    else:raise ValueError("Not an approved source exposure sample weighting")
    terms=["site_light"]+(["doy_30"] if include_doy else [])
    if yearFE=="species_year":
        second=d.species_year_key.to_numpy()
    elif yearFE=="year":
        second=d.year.to_numpy()
    else:raise ValueError("Wrong source time FE")
    yy=d[y].to_numpy(float)
    xx=d[terms].to_numpy(float)
    z,fe=weighted_panel_partial_out(
        np.column_stack((yy,xx)),
        (d.panel_key.to_numpy(),second),w)
    ytilde=z[:,0];xtilde=z[:,1:]
    p=len(terms)
    mat=(xtilde*w[:,None]).T@xtilde
    condition=float(np.linalg.cond(mat))
    if not np.isfinite(condition) or condition>1e12 or np.linalg.matrix_rank(mat)<p:
        return {"status":"source_light_or_doy_unidentified_after_original_FEs",
                "condition_number":condition,
                "source_true_plots":G,
                "source_rows":N,"fe":fe}
    bread=np.linalg.inv(mat)
    betas=bread@((xtilde*w[:,None]).T@ytilde)
    u=ytilde-xtilde@betas
    cid=pd.factorize(d.location_id,sort=True)[0]
    scores=np.zeros((G,p),dtype=float)
    np.add.at(scores,cid,xtilde*(w*u)[:,None])
    vcov=bread@(scores.T@scores)@bread
    vcov*=G/(G-1)*N/(N-p)
    ses=np.sqrt(np.clip(np.diag(vcov),0,None))
    require(np.isfinite(ses).all() and np.any(ses>0),
            "No genuine plot-level cluster uncertainty identified")
    q=float(t.ppf(.975,G-1))
    coeff={}
    for i,label in enumerate(terms):
        coeff[label]={"beta":float(betas[i]),
                      "genuine_plot_cluster_SE":float(ses[i]),
                      "nominal_95CI":[float(betas[i]-q*ses[i]),float(betas[i]+q*ses[i])],
                      "nominal_unadjusted_p":float(2*t.sf(abs(betas[i]/ses[i]),G-1))
                         if ses[i]>0 else None}
    rsq=1-float(np.sum(w*u*u)/max(np.sum(w*ytilde*ytilde),1e-12))
    return {
       "status":"estimated_exploratory_within_positive_recorded_plant_cover",
       "recorded_outcome":y,
       "source_plot_species_surveyID_season_FE":True,
       "source_original_year_FE":yearFE,
       "cluster_real_physical_NPMS_plot":True,
       "source_weighting":equal,
       "source_covariates":terms,
       "actual_original_recorded_cover_rows":N,
       "actual_independent_source_plots":G,
       "source_distinct_taxa":int(d.species.nunique()),
       "source_distinct_original_years":int(d.year.nunique()),
       "source_original_year_range":[int(d.year.min()),int(d.year.max())],
       "source_panel_FE_support":fe,
       "condition_number_source_covariates":condition,
       "source_within_FE_weighted_R2_descriptive":rsq,
       "terms":coeff,
       "positive_only_missing_species_not_absence":True,
       "not_causal_ALAN_on_plant_cover_or_actual_nocturnal_pollination":True}

def run(source,receipt,glopl,out):
    df,audit=read_original(source,receipt)
    g=pd.read_csv(glopl,encoding="latin-1",low_memory=False,dtype=str)
    require(len(g)==2969 and sha(glopl)==
            "264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de",
            "Real original GloPL independent study source not verified")
    name=g.Species_accepted_names.fillna("").str.replace("_"," ",regex=False).str.strip()
    g["species"]=name
    g["study_key"]=np.where(g.DOI.fillna("").str.strip().ne(""),
              g.DOI.fillna("").str.casefold(),g.Author.fillna("")+"|"+g.Year.fillna(""))
    studies=g.groupby("species").study_key.nunique()
    glopl=set(studies.index)
    reliable=set(studies.index[studies>=2])
    df["source_GloPL_recorded"]=df.species.isin(glopl)
    df["source_GloPL_two_studies"]=df.species.isin(reliable)
    base,selected=choose(df)
    variants={
       "primary_original_cover_percentage_pp_source_plot_equal_with_species_year_FE":
          fit_source(base,"recorded_cover_pp","plot_equal","species_year",True),
       "original_cover_log1p_percentage":
          fit_source(base,"recorded_cover_log","plot_equal","species_year",True),
       "source_record_year_equal":
          fit_source(base,"recorded_cover_pp","record_equal","species_year",True),
       "source_equal_plot_species_panel":
          fit_source(base,"recorded_cover_pp","plot_species_equal","species_year",True),
       "year_FE_no_species_specific_year_control":
          fit_source(base,"recorded_cover_pp","plot_equal","year",True),
       "source_actual_DOY_not_adjusted":
          fit_source(base,"recorded_cover_pp","plot_equal","species_year",False),
    }
    bands={}
    for threshold in (.1,.3):
        z,diag=choose(df,source_min_range=threshold)
        bands[str(threshold)]={"source_eligibility":diag,
                          "model":fit_source(z,"recorded_cover_pp","plot_equal",
                                             "species_year",True)}
    method={}
    for survey in ("87","154","155"):
        z,diag=choose(df,only_survey_id=survey)
        method[survey]={"source_eligibility":diag,
                         "model":fit_source(z,"recorded_cover_pp","plot_equal",
                                             "species_year",True)}
    z,doydiag=choose(df,max_date_deviation=30)
    yearstatus,statdiag=choose(df,source_status_stable=True)
    sourcePL={}
    for tag,field in (
        ("exact_GloPL_original_supplementation_species","source_GloPL_recorded"),
        ("two_or_more_independently_studied_GloPL_species","source_GloPL_two_studies")
    ):
        z=base.loc[base[field]].copy()
        z,newdiag=choose(z)
        sourcePL[tag]={"source_eligibility":newdiag,
                       "study_supported_species":int(z.species.nunique()),
                       "actual_source_physical_plots":int(z.location_id.nunique()),
                       "model":fit_source(z,"recorded_cover_pp","plot_equal",
                                           "species_year",True)}
    # A methodological negative-control response: number of *positive source*
    # cover records at this plant-site-year is an outcome ascertainment proxy,
    # not total recording effort or plot-wide species richness.
    variants["positive_plant_recording_opportunity_negative_control"]=fit_source(
          base,"log_recorded_visits","plot_equal","species_year",True)
    variants["same_visit_DOY_plusminus30days"]=fit_source(
          z if False else choose(df,max_date_deviation=30)[0],
          "recorded_cover_pp","plot_equal","species_year",True)
    variants["same_source_record_status_across_years"]=fit_source(
          yearstatus,"recorded_cover_pp","plot_equal","species_year",True)

    output={
       "schema":SCHEMA,
       "source_original_NPMS_DOI":"10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8",
       "source_verified_genuine_GloPL_DOI":"10.5061/dryad.dt437",
       "source_corrected_2014_2024_previous_PCNL_DOI":"10.5281/zenodo.17013414",
       "source_original_positive_record_cover_CSV_SHA256":sha(source),
       "original_all_positive_source_plant_record_rowyears":len(df),
       "source_primary_eligibility":selected,
       "true_source_NPMS_eligible_study_years":[EARLIEST,LATEST],
       "original_DOMIN_midpoints_are_fraction_0_to_1":True,
       "source_original_plant_species_absence_not_imputed_from_nonrecord":True,
       "source_repeated_plots_not_random_assignment_of_light":True,
       "source_preyear_light_not_local_LUX_or_flower_light":True,
       "source_cover_response_measured_only_conditional_on_NPMS_plant_record":True,
       "source_GloPL_independent_study_PL_not_NPMS_site_measured_PL":True,
       "source_covariate_original_day_of_year_recorded":True,
       "original_five_model_specifications_and_sampling_robustness":variants,
       "source_light_range_sensitivity":bands,
       "source_same_original_surveyID_season_strata":method,
       "external_GloPL_study_replication_scope":sourcePL,
       "source_phenology_matched_Doy_30day_eligibility":doydiag,
       "source_record_status_constant_across_year_eligibility":statdiag,
       "unmeasured_habitat_management_weather_observer_switching_confounding":True,
       "not_measured_real_moth_pollination_or_causal_plant_decline":True,
       "all_p_values_nominal_unadjusted_postexploratory":True}
    out=Path(out)
    require(not out.exists(),"No overwrite original model")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    keep=variants[
      "primary_original_cover_percentage_pp_source_plot_equal_with_species_year_FE"]
    print("ALAN_TRUE_NPMS_SOURCE_PLOT_POSITIVE_COVER_VS_PREALAN",
          json.dumps({"primary_source_eligibility":selected,
                      "primary":keep,
                      "year_only":variants["year_FE_no_species_specific_year_control"],
                      "GloPL_two_studies":sourcePL[
                         "two_or_more_independently_studied_GloPL_species"]},
                     sort_keys=True),flush=True)
    return output

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("source","receipt","glopl","out"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.receipt,a.glopl,a.out)

