#!/usr/bin/env python3
"""Original 703-species ALAN × directly sampled same-year WorldPop population.

Exposure: prior-to-photo-year Li-Cao PCNL mean log1p, observational.
Outcome: species photographically recorded or unrecorded in original plant
photo opportunity 0.05-degree cells. Nonrecord != biological absence.

A single matched source cohort is used for both baseline (light+effort)
and conditional (light+effort+log1p(pop density)) models. Historical R2024B
UC population was measured in exactly the same WorldPop 30 arcsecond pixel,
not an arbitrary other plant photo inside the same 0.05-degree grid.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, load_sites, load_labels, original_accepted_species, require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import (
    direct_plant_traits
)
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import (
    prepare_support, island_fold, percentile
)
from scripts.analysis.alan_mois_worldpop_exact_pixel_support_20261115 import join

N_BOOT=400
SEED=20261115

def within_year_island_photo_design(pop_sites):
    p=pop_sites.loc[pop_sites.pop_density.notna()].copy()
    require(not p.duplicated(KEY).any(),"Population subset source key duplicates")
    p["logpop"]=np.log1p(p.pop_density.to_numpy(float))
    p["loglight"]=p[LIGHT].astype(float)
    p["effort"]=p.log_obs.astype(float)
    g=p.groupby(["island","year"],observed=True).loglight.agg(
        size="size",min="min",max="max",median="median").reset_index()
    g=g.loc[g["size"].ge(6)&(g["max"]-g["min"]).ge(.02)]
    p=p.merge(g[["island","year","median"]],
              on=["island","year"],how="inner",validate="many_to_one")
    p["hi"]=(p.loglight>p["median"]).astype("int8")
    halves=p.groupby(["island","year","hi"],observed=True).size().unstack(
        "hi",fill_value=0)
    for k in (0,1):
        if k not in halves:halves[k]=0
    ggood=halves.loc[halves[0].ge(2)&halves[1].ge(2)].reset_index()[
        ["island","year"]]
    p=p.merge(ggood,on=["island","year"],how="inner",validate="many_to_one")
    p=p.reset_index(drop=True)
    p["sid"]=np.arange(len(p),dtype=np.int32)
    p["iy"],_=pd.factorize(p.island.astype(str).str.cat(
       p.year.astype(str),sep="|"),sort=True)
    p["fold"]=p.island.map(island_fold)
    require(len(p)>0,"No valid independent matched population gradients")
    # Cohort-wide within-island/year partial correlations; site-level controls only.
    covs=p[["loglight","logpop","effort"]]
    demeaned=covs-p.groupby("iy",observed=True)[list(covs)].transform("mean")
    corr=demeaned.corr()
    return p,{
      "sameyear_pixel_matched_original_site_years_before_gradient":int(pop_sites.pop_density.notna().sum()),
      "sameyear_pixel_matched_source_islands_before_gradient":int(
         pop_sites.loc[pop_sites.pop_density.notna(),"island"].nunique()),
      "eligible_photo_site_years_after_pop_and_light_gradient":len(p),
      "eligible_islands_after_pop_and_light_gradient":int(p.island.nunique()),
      "eligible_island_by_year_strata":int(p.iy.nunique()),
      "within_island_year_logPCNL_logpop_correlation":float(corr.loc["loglight","logpop"]),
      "within_island_year_logPCNL_effort_correlation":float(corr.loc["loglight","effort"]),
      "within_island_year_logpop_effort_correlation":float(corr.loc["logpop","effort"]),
      "historical_NTL_cannot_be_called_local_light_treatment":True,
    }

def candidate_species_photo_support(p,flora,autonomy):
    rows=flora[KEY+["species"]].merge(
        p[KEY+["sid","iy","hi"]],
        on=KEY,how="inner",validate="many_to_one")
    rows=rows.loc[rows.species.isin(autonomy.index)].copy()
    strata_n=p.groupby("iy",observed=True).size().rename("n_sites")
    st=(rows.groupby(["species","iy"],observed=True)
         .agg(n_positive=("sid","size"),n_low=("hi",lambda x:int(x.eq(0).sum())),
              n_high=("hi","sum")).reset_index())
    st["n_sites"]=st.iy.map(strata_n)
    st=st.loc[st.n_positive.lt(st.n_sites)].copy()
    site_ids={sp:g.sid.to_numpy(np.int64) for sp,g in
              rows.groupby("species",observed=True)}
    strata={sp:g.iy.to_numpy(int) for sp,g in
            st.groupby("species",observed=True)}
    return site_ids,strata,{
        "population_eligible_original_species_positive_rows":len(rows),
        "documented_candidates_with_at_least_one_photo":int(rows.species.nunique()),
        "documented_candidates_with_at_least_one_informative_island_year":len(strata)
    }

def fit_species(p,positive,iy_set,min_positive=12,min_strata=4,min_islands=2):
    s=p.loc[p.iy.isin(iy_set)].copy()
    if len(s)<30:return None,"TOO_FEW_PHOTO_SITE_OPPORTUNITIES"
    y=np.isin(s.sid.to_numpy(np.int64),positive,assume_unique=True).astype(float)
    npos=int(y.sum())
    nisland=int(s.island.nunique())
    nstrata=int(s.iy.nunique())
    if npos<min_positive or nstrata<min_strata or nisland<min_islands:
        return None,"INSUFFICIENT_SPECIES_POP_MATCHED_SUPPORT"
    if npos>=len(s):return None,"NO_RECORDED_VS_NONRECORDED_CONTRAST"
    cols=["loglight","effort","logpop"]
    cols_y=s[cols].copy()
    cols_y["y"]=y
    mu=cols_y.groupby(s.iy).transform("mean")
    x=(cols_y-mu).to_numpy(float)
    yy=x[:,-1]
    X=x[:,:3]
    # Each island contributes unit total source sampling weight for that species.
    codes,levels=pd.factorize(s.island,sort=True)
    w=1/np.bincount(codes)[codes].astype(float)
    scales=np.sqrt(np.sum(w[:,None]*X**2,axis=0)/np.sum(w))
    if not np.isfinite(scales).all() or np.min(scales)<1e-9:
        return None,"NO_COVARIATE_RESIDUAL_VARIATION"
    z=X/scales
    fit={}
    for label,ncov in (("same_matched_cohort_light_effort",2),
                       ("conditional_light_effort_population",3)):
        XX=z[:,:ncov]
        M=(XX.T*w)@XX
        if np.linalg.matrix_rank(M)<ncov or np.linalg.cond(M)>1e8:
            return None,"LIGHT_POP_EFFORT_COLLINEAR"
        bread=np.linalg.inv(M)
        beta=bread@(XX.T@(w*yy))
        raw=beta/scales[:ncov]
        if not np.isfinite(raw).all():return None,"NONFINITE_COVARIATE_BETA"
        fit[label]=raw
    result={
      "n_positive":npos,"n_source_sites":len(s),
      "n_islands":nisland,"n_strata":nstrata,
      "beta_L_without_population_pp":float(100*fit["same_matched_cohort_light_effort"][0]),
      "beta_L_conditional_population_pp":float(100*fit["conditional_light_effort_population"][0]),
      "beta_population_conditional_light_pp":float(100*fit["conditional_light_effort_population"][2]),
      "beta_effort_conditional_population_pp":float(100*fit["conditional_light_effort_population"][1]),
      "population_changes_ALAN_beta_pp":float(100*(
            fit["conditional_light_effort_population"][0]-
            fit["same_matched_cohort_light_effort"][0])),
    }
    # Independent island-cluster source partial-light SE, conditional model.
    if nisland>=10:
        XX=z
        bb=fit["conditional_light_effort_population"]*scales
        resid=yy-XX@bb
        bread=np.linalg.inv((XX.T*w)@XX)
        score=np.zeros((nisland,3),float)
        np.add.at(score,codes,XX*(w*resid)[:,None])
        dof=len(s)-nstrata-3
        if dof>0:
            vc=((nisland/(nisland-1))*((len(s)-1)/dof))*bread@(score.T@score)@bread
            se=np.sqrt(np.maximum(np.diag(vc),0))/scales
            cutoff=student_t.ppf(.975,nisland-1)
            result["se_L_popadjusted_pp"]=float(100*se[0])
            result["ci_L_popadjusted_pp"]=[
                float(100*(fit["conditional_light_effort_population"][0]-cutoff*se[0])),
                float(100*(fit["conditional_light_effort_population"][0]+cutoff*se[0]))]
            result["se_population_conditional_L_pp"]=float(100*se[2])
    # Island-specific sufficient statistics for shared-block resampling; store
    # in RAM only, never export species×island data or trait crosswalk.
    Xraw=X
    stats_A=np.zeros((nisland,3,3),float)
    stats_B=np.zeros((nisland,3),float)
    for k in range(nisland):
        idx=np.flatnonzero(codes==k)
        a=Xraw[idx]
        ww=w[idx]
        stats_A[k]=(a.T*ww)@a
        stats_B[k]=a.T@(ww*yy[idx])
    # Return result and IDs/sufficient stats for geographic bootstrap.
    return (result,levels.astype(str).to_numpy(),
            stats_A,stats_B),None

def aggregate(fits):
    d=pd.DataFrame([{k:v for k,v in row.items() if k not in ("species","genus")}
                    for row in fits])
    out={"fit_species":len(d)}
    for level,part in (("all",d),("autonomous",d.loc[d.autonomy.eq("autonomous")]),
                       ("absent",d.loc[d.autonomy.eq("absent")])):
        out[level]={
           "n_species":len(part),
           "island_support_ge10":int(part.n_islands.ge(10).sum()),
           "island_support_ge30":int(part.n_islands.ge(30).sum()),
           "beta_L_without_pop_pp":percentile(part.beta_L_without_population_pp),
           "beta_L_given_pop_pp":percentile(part.beta_L_conditional_population_pp),
           "beta_pop_given_L_pp":percentile(part.beta_population_conditional_light_pp),
           "delta_beta_L_after_pop_pp":percentile(part.population_changes_ALAN_beta_pp),
        }
    for key in ("beta_L_without_population_pp",
                "beta_L_conditional_population_pp",
                "beta_population_conditional_light_pp"):
        out[key+"_autonomous_minus_absent"]=float(
            d.loc[d.autonomy.eq("autonomous"),key].mean()-
            d.loc[d.autonomy.eq("absent"),key].mean())
    return out

def block_bootstrap(fits,universe,nrep=N_BOOT):
    # Cohort restricted ex ante to species with >=10 independent islands and
    # sufficient original matched covariate support; no p-value selection.
    eligible=[r for r in fits if r["fit"]["n_islands"]>=10]
    yes=sum(x["autonomy"]=="autonomous" for x in eligible)
    no=len(eligible)-yes
    if yes<25 or no<25:
        return {"status":"NOT_ESTIMABLE_25_SPECIES_PER_SELFING_GROUP",
                "species":len(eligible),"yes":yes,"no":no}
    universe=list(universe)
    lookup={name:i for i,name in enumerate(universe)}
    n=len(eligible);g=len(universe)
    mats=np.zeros((n,g,3,3),float)
    vecs=np.zeros((n,g,3),float)
    aut=np.zeros(n,bool)
    for i,t in enumerate(eligible):
        aut[i]=t["autonomy"]=="autonomous"
        ind=np.asarray([lookup[j] for j in t["islands"]],dtype=int)
        mats[i,ind,:,:]=t["A"]
        vecs[i,ind,:]=t["B"]
    # Baseline full-clone verified: this is the same 3-regressor original
    # within-photo-stratum WLS, with exact same per-island equal weighting.
    aa=mats.sum(axis=1)
    bb=vecs.sum(axis=1)
    recovered=np.linalg.solve(aa,bb)
    target=np.array([t["fit"]["beta_L_conditional_population_pp"] for t in eligible])
    require(np.max(abs(100*recovered[:,0]-target))<1e-7,
            "Island block sufficient stats do not reconstruct source fit")
    rng=np.random.default_rng(SEED)
    boot=np.full((nrep,2),np.nan)
    for b in range(nrep):
        mult=np.bincount(rng.integers(0,g,g),minlength=g).astype(float)
        aa=np.einsum("g,ngjk->njk",mult,mats,optimize=True)
        bb=np.einsum("g,ngj->nj",mult,vecs,optimize=True)
        try:
            beta=np.linalg.solve(aa,bb)
        except np.linalg.LinAlgError:
            continue
        beta=beta*100
        if not np.isfinite(beta).all():continue
        boot[b,0]=beta[aut,0].mean()-beta[~aut,0].mean()
        boot[b,1]=beta[aut,2].mean()-beta[~aut,2].mean()
    vals=boot[np.isfinite(boot).all(axis=1)]
    if len(vals)<nrep*.8:
        return {"status":"NOT_ESTIMABLE_BOOTSTRAP_RANK_STABILITY",
                "n_species":n,"n_valid_draws":len(vals)}
    return {
      "status":"CONDITIONAL_PHOTO_OBSERVATION_ISLAND_BLOCK_BOOTSTRAP",
      "n_species_ge10_islands":n,"autonomous":int(aut.sum()),
      "absent":int((~aut).sum()),"n_source_islands":g,
      "n_draws_requested":nrep,"n_valid":len(vals),
      "conditional_ALAN_slope_group_diff_pp":float(
            target[aut].mean()-target[~aut].mean()),
      "conditional_population_group_diff_pp":float(
            100*recovered[aut,2].mean()-100*recovered[~aut,2].mean()),
      "ALAN_group_diff_95percentile_ci_pp":[float(x) for x in
             np.quantile(vals[:,0],[.025,.975])],
      "population_group_diff_95percentile_ci_pp":[float(x) for x in
             np.quantile(vals[:,1],[.025,.975])],
      "global_island_resampling_conditional_on_nonrandom_population_coverage":True
    }

def run(site_file,globi,groups,ledger,pop_sources,out):
    require(not out.exists(),"Refuse overwrite")
    sites=load_sites(site_file)
    flora=original_accepted_species(groups,sites,load_labels(globi))
    trait=direct_plant_traits(ledger,high_only=False)
    names=pd.DataFrame({"species":flora.species.drop_duplicates()})
    names["autonomy"]=names.species.map(trait.autonomous_selfing_capacity)
    # Exact original eligible 703 before any second-source population selection.
    _,_,_,original703,support=prepare_support(sites,flora,names,span=.02)
    require((len(original703),support["eligible_islands"])==(703,243),
            "Original 703 species denominator drift")
    original_aut=names.set_index("species").autonomy.loc[sorted(original703)]
    match,popframe=join(sites,pop_sources,return_data=True)
    p,design=within_year_island_photo_design(popframe)
    site_ids,strata,counts=candidate_species_photo_support(
        p,flora,original_aut)
    fits=[];raw=[];failed={}
    for j,sp in enumerate(sorted(original703)):
        if sp not in strata or sp not in site_ids:
            failed["NO_INFORMATIVE_SOURCE_STRATA"]=failed.get("NO_INFORMATIVE_SOURCE_STRATA",0)+1
            continue
        ans,reason=fit_species(p,site_ids[sp],strata[sp])
        if ans is None:
            failed[reason]=failed.get(reason,0)+1
            continue
        fit,island_ids,A,B=ans
        row={"species":sp,"autonomy":str(original_aut[sp]),**fit}
        fits.append(row)
        raw.append({"fit":fit,"autonomy":str(original_aut[sp]),
                    "islands":island_ids,"A":A,"B":B})
        if (j+1)%100==0: print("ALAN_POPULATION_SOURCE_SPECIES_PROGRESS",j+1,flush=True)
    require(len(fits)>0,"No population-adjusted species actually estimable")
    out_summary={
      "schema":"ALAN_703_selfing_species_exact_annual_WorldPop_population_conditional_photo_slope_v1",
      "status":"ORIGINAL_SOURCE_SUBCOHORT_OBSERVATIONAL_POPULATION_CONTROL_NOT_CAUSAL",
      "sources":{"original_photo_site_years":77647,"original_islands":562,
          "original_source_species_site_years":2012925,
          "fixed_original703_source_trait_taxa":703,
          "MoIS_annual_WorldPop_G2_sameyear_samepixel_R2024B_UC":True,
          "no_new_population_imputation":True},
      "population_source_join":match,
      "matched_light_design":design,
      "species_support":{**counts,"original_fixed_703":len(original703),
                         "fit_with_population":len(fits),
                         "fitted_selfing_yes":int(sum(x["autonomy"]=="autonomous" for x in fits)),
                         "fitted_selfing_no":int(sum(x["autonomy"]=="absent" for x in fits)),
                         "not_fitted_source_reasons":failed},
      "conditional_model_comparison":aggregate(fits),
      "independent_island_bootstrap":block_bootstrap(raw,p.island.unique()),
      "limitations":{
        "population_observed_at_only_30pct_of_original_photo_sites":True,
        "MOIS_photo_cohort_selection_bias_unresolved":True,
        "population_and_NTL_same_original_calendar_observation_year_but_NTL_is_prior_mean":True,
        "no_climate_or_landuse_adjustment":True,
        "photo_nonrecord_not_true_absence":True,
        "no_field_illuminance_or_light_experiment":True,
        "no_true_pollinator_or_seed_reproduction_data":True,
        "no_species_level_private_results_exported":True,
        "multiple_exploratory_model_versions_not_confirmation":True
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(out_summary,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_WORLDPOP_CONDITIONAL_703_SOURCE",json.dumps({
      "population_join_count":match["matched_usable_site_years"],
      "matched_light_design":design,
      "species_support":out_summary["species_support"],
      "conditional_model":out_summary["conditional_model_comparison"],
      "bootstrap":out_summary["independent_island_bootstrap"],
    },ensure_ascii=False,allow_nan=False),flush=True)
    return out_summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for opt in ("sites","globi","groups","ledger","pop_sources","out"):
        p.add_argument("--"+opt,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.pop_sources,a.out)
