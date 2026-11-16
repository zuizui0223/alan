#!/usr/bin/env python3
"""Independent JRC GHSL urban built-up control for ALAN plant photo slopes.

Pre-registered two sources:
- GHS-BUILT-S epoch 2015, source prior to all photographed plants: PRIMARY.
- Epoch 2020: sensitivity only (later than observations in 2019).

This analysis estimates photographically-recorded plant species rather than
measured plant abundance/occupancy, seed or actual pollination. Even the
GHSL-adjusted result cannot identify causal lamp light effects.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites,load_labels,original_accepted_species,require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import prepare_support,percentile
from scripts.analysis.alan_species_worldpop_conditional_20261115 import (
    within_year_island_photo_design,candidate_species_photo_support
)
from scripts.analysis.alan_species_worldclim_population_20261115 import matching_climate_gradient,fit_five
from scripts.analysis.alan_species_full_worldpop_climate_20261116 import source_annual_population
from scripts.analysis.alan_worldpop_full_77647_sites_20261116 import file_sha256

N_B=400
SOURCE_GHSL_NAME="alan_77647_original_ghsl_build_2015_2020_source_sitevalues.csv.gz"
ORDER=("loglight","effort","logpop","bio1","logprec")
BASELEN=len(ORDER)

def authentic_ghsl_source(source,climate):
    source=Path(source)
    rec=json.loads((source/"source_receipt.json").read_text())
    f=source/SOURCE_GHSL_NAME
    require(rec["schema"]=="ALAN_77k_photo_original_GHSL_BUILT_S_R2023A_30asec_epochs2015_2020_v1"
            and rec["original_ALAN_photo_site_years"]==77647
            and rec["source_2015_independent_prior_urban_builtup"]
            and rec["source_2020_sensitivity_only"]
            and rec["no_builtup_source_nodata_imputed_zero"],
            "Wrong original GHSL raster sample source contract")
    require(file_sha256(f)==rec["output_site_values_sha256"],
            "GHSL site values SHA wrong")
    raw=pd.read_csv(f,dtype={"site_id":np.int32})
    require(raw.site_id.is_unique and len(raw)==77647,
            "GHSL original ALAN source site identities must be one-to-one")
    z=climate.merge(raw[["site_id",
        "ghsl_built_2015_fraction","ghsl_built_2020_fraction"]],
        on="site_id",how="left",validate="one_to_one")
    require(len(z)==len(climate),"GHSL source sample join expanded/dropped original")
    for e in (2015,2020):
        z[f"logbuilt{e}"]=np.log1p(100*np.maximum(
            0,z[f"ghsl_built_{e}_fraction"].to_numpy(float)))
    both=z.logbuilt2015.notna()&z.logbuilt2020.notna()
    z=z.loc[both].copy()
    require(len(z)>1000 and z.island.nunique()>30,
            "Insufficient GHSL source complete-case coverage")
    # Reapply the SAME original within-island/year exposure and presence
    # opportunity screen after any genuine raster NODATA exclusion.
    g=z.groupby(["island","year"],observed=True).loglight.agg(
        n="size",lo="min",hi="max",med="median")
    ok=g.loc[g.n.ge(6)&(g.hi-g.lo).ge(.02)]
    z=z.merge(ok[["med"]].reset_index(),on=["island","year"],
              how="inner",validate="many_to_one")
    z["higher"]=z.loglight.gt(z.med)
    groups=z.groupby(["island","year","higher"],observed=True).size().unstack(
        "higher",fill_value=0)
    for b in (False,True):
        if b not in groups:groups[b]=0
    matched=groups.loc[groups[False].ge(2)&groups[True].ge(2)].reset_index()[
        ["island","year"]]
    z=z.merge(matched,on=["island","year"],validate="many_to_one").reset_index(drop=True)
    z["sid"]=np.arange(len(z),dtype=np.int32)
    z["iy"],_=pd.factorize(z.island.astype(str).str.cat(
        z.year.astype(str),sep="|"),sort=True)
    cov=list(ORDER)+["logbuilt2015","logbuilt2020"]
    residual=z[cov]-z.groupby("iy",observed=True)[cov].transform("mean")
    sd=residual.std()
    return z,{
        "source_ghsl_2015_raster_sha256":
            rec["source"]["2015"]["sha256_original_unmodified_tif"],
        "source_ghsl_2020_raster_sha256":
            rec["source"]["2020"]["sha256_original_unmodified_tif"],
        "original_population_climate_source_siteyears_before_build_control":len(climate),
        "source_2015_and_2020_builtup_numeric_complete":int(both.sum()),
        "eligible_same_site_years_after_builtin_source_and_light_support":len(z),
        "source_islands":int(z.island.nunique()),
        "source_island_years":int(z.iy.nunique()),
        "original_prior_builtup_2015_source":True,
        "epoch2020_future_exposure_for_observations2019_sensitivity_only":True,
        "median_built2015_fraction_source":
           float(z.ghsl_built_2015_fraction.median()),
        "median_built2020_fraction_source":
           float(z.ghsl_built_2020_fraction.median()),
        "built2015_zero_source_siteyears":int(z.ghsl_built_2015_fraction.eq(0).sum()),
        "built2020_zero_source_siteyears":int(z.ghsl_built_2020_fraction.eq(0).sum()),
        "global_within_islandyear_residual_corr":{
           k:{p:float(residual[k].corr(residual[p])) for p in cov}
              for k in cov},
        "global_within_islandyear_residual_sd":{
           k:float(sd[k]) for k in cov},
        "original_landcover_not_plant_natural_habitat":True
    }

def fit_source_species(sites,ids,informative_iy):
    s=sites.loc[sites.iy.isin(informative_iy)].copy()
    if len(s)<30:return None,"INSUFFICIENT_SOURCE_SITE_YEARS"
    y=np.isin(s.sid.to_numpy(np.int64),ids,assume_unique=True).astype(float)
    G=int(s.island.nunique());K=int(s.iy.nunique());npos=int(y.sum())
    if npos<12 or K<4 or G<2:
        return None,"INSUFFICIENT_SOURCE_POSITIVE_ISLAND_SUPPORT"
    if npos>=len(s):return None,"MISSING_SOURCE_PHOTO_NONRECORD"
    fullcols=list(ORDER)+["logbuilt2015","logbuilt2020"]
    original=s[fullcols].copy()
    original["y"]=y
    within=(original-original.groupby(s.iy).transform("mean")).to_numpy(float)
    X=within[:,:7];yy=within[:,-1]
    codes,islands=pd.factorize(s.island,sort=True)
    w=1/np.bincount(codes)[codes].astype(float)
    norms=np.sqrt((w[:,None]*X**2).sum(axis=0)/w.sum())
    if np.min(norms)<1e-9 or not np.isfinite(norms).all():
        return None,"NO_SOURCE_WITHIN_ISLANDYEAR_BUILDING_OR_POPULATION_VARIATION"
    Z=X/norms
    models={}
    slices={"baseline5":[0,1,2,3,4],
            "ghsl2015":[0,1,2,3,4,5],
            "ghsl2020":[0,1,2,3,4,6]}
    for label,columns in slices.items():
        XX=Z[:,columns]
        M=(XX.T*w)@XX
        if np.linalg.matrix_rank(M)!=len(columns) or np.linalg.cond(M)>1e8:
            return None,"RANK_OR_COLLINEARITY_GHSL_POP_LIGHT_CLIMATE"
        sol=np.linalg.solve(M,XX.T@(w*yy))
        raw=sol/norms[columns]
        if not np.isfinite(raw).all():
            return None,"NONFINITE_GHSL_SOURCE_SLOPE"
        models[label]=raw
    q={
      "n_islands":G,"n_island_year_strata":K,
      "n_photo_species_positive":npos,
      "n_original_site_years":len(s),
      "L_baseline_pp":float(100*models["baseline5"][0]),
      "L_built2015_pp":float(100*models["ghsl2015"][0]),
      "L_built2020_pp":float(100*models["ghsl2020"][0]),
      "built2015_pp_per_log_percent_built":float(100*models["ghsl2015"][5]),
      "built2020_pp_per_log_percent_built":float(100*models["ghsl2020"][5]),
      "L_delta_post_build2015_pp":float(
          100*(models["ghsl2015"][0]-models["baseline5"][0])),
    }
    # Method comparison, exact original 5-variable fit on identical sites.
    prior,why=fit_five(s,ids,informative_iy)
    require(why is None,"Original population+climate five-variable source failed")
    v,*_=prior
    require(abs(q["L_baseline_pp"]-v["beta_L_with_climate_pp"])<1e-7,
            "Same original population+climate beta not reproduced")
    # Per-island score blocks for actual island re-sampling and ALL coefficient
    # refits. DO NOT export these original species×island matrices.
    A=np.zeros((G,7,7),dtype=float)
    B=np.zeros((G,7),dtype=float)
    for i in range(G):
        ix=np.flatnonzero(codes==i)
        x=X[ix];ww=w[ix]
        A[i]=(x.T*ww)@x
        B[i]=x.T@(ww*yy[ix])
    return (q,islands.astype(str).to_numpy(),A,B),None

def groups_summary(fits):
    d=pd.DataFrame(fits)
    out={"n_model_fitted":len(d)}
    cols=("L_baseline_pp","L_built2015_pp","L_built2020_pp",
          "built2015_pp_per_log_percent_built","built2020_pp_per_log_percent_built",
          "L_delta_post_build2015_pp")
    for name,mask in (
       ("all",np.ones(len(d),dtype=bool)),
       ("selfing_autonomous",d.autonomy.eq("autonomous").to_numpy(bool)),
       ("selfing_absent",d.autonomy.eq("absent").to_numpy(bool))):
        z=d.loc[mask].copy()
        out[name]={
          "n_species":len(z),
          "n_with_ge10_independent_islands":int(z.n_islands.ge(10).sum()),
          "n_with_ge20_independent_islands":int(z.n_islands.ge(20).sum()),
          "n_with_ge30_independent_islands":int(z.n_islands.ge(30).sum()),
          **{c:percentile(z[c]) for c in cols}
        }
    for c in cols:
        out[c+"_autonomous_minus_absent"]=float(
            d.loc[d.autonomy.eq("autonomous"),c].mean()-
            d.loc[d.autonomy.eq("absent"),c].mean())
    return out

def island_bootstrap(all_models,source_islands,epoch="ghsl2015"):
    label="ghsl2015" if epoch=="ghsl2015" else "ghsl2020"
    chosen=[x for x in all_models if x["q"]["n_islands"]>=20]
    yes=sum(x["autonomy"]=="autonomous" for x in chosen)
    no=len(chosen)-yes
    if yes<25 or no<25:
        return {"status":"INSUFFICIENT_20ISLAND_GROUP_SUPPORT",
                "n_species":len(chosen),"autonomous":yes,"absent":no}
    universe=list(source_islands)
    ix={str(x):i for i,x in enumerate(universe)}
    n=len(chosen);g=len(universe)
    mats=np.zeros((n,g,7,7),float)
    vecs=np.zeros((n,g,7),float)
    aut=np.zeros(n,bool)
    for i,row in enumerate(chosen):
        aut[i]=row["autonomy"]=="autonomous"
        ind=np.asarray([ix[x] for x in row["islands"]],dtype=int)
        mats[i,ind,:,:]=row["A"]
        vecs[i,ind,:]=row["B"]
    cols=np.array([0,1,2,3,4,5 if label=="ghsl2015" else 6],int)
    mats=mats[:,:,cols,:][:,:,:,cols]
    vecs=vecs[:,:,cols]
    def solve(weight):
        m=np.einsum("g,ngjk->njk",weight,mats,optimize=True)
        b=np.einsum("g,ngj->nj",weight,vecs,optimize=True)
        return np.linalg.solve(m,b[...,None])[...,0]
    first=solve(np.ones(g))
    original=np.array([x["q"]["L_built2015_pp" if label=="ghsl2015" else "L_built2020_pp"] for x in chosen])
    require(np.isfinite(first).all() and
           abs(100*first[:,0]-original).max()<1e-7,
           "Independent original GHSL source sufficient stats fail to reproduce first stage")
    rng=np.random.default_rng(20261116)
    boot=np.full((400,2),np.nan)
    for z in range(400):
        counts=np.bincount(rng.integers(0,g,g),minlength=g).astype(float)
        try:
            b=solve(counts)*100
        except np.linalg.LinAlgError:continue
        if not np.isfinite(b).all():continue
        boot[z,0]=b[aut,0].mean()-b[~aut,0].mean()
        boot[z,1]=b[aut,5].mean()-b[~aut,5].mean()
    valid=boot[np.isfinite(boot).all(axis=1)]
    if len(valid)<320:
        return {"status":"NOT_ESTIMABLE_ISLAND_BLOCK_RANK",
                "n_valid":len(valid),"n_requested":400,
                "n_species":n,"n_independent_islands":g}
    return {
      "status":"ORIGINAL_SOURCE_2015_PREPHOTO_BUILTUP_ISLAND_BOOTSTRAP"
         if epoch=="ghsl2015" else "EPOCH_2020_BUILTUP_SENSITIVITY_BOOTSTRAP",
      "n_species_with_ge20_islands":n,
      "autonomous":yes,"absent":no,
      "n_sampling_islands":g,
      "source_epoch":2015 if epoch=="ghsl2015" else 2020,
      "n_valid":len(valid),"n_requested":400,
      "L_autonomous_minus_absent_pp":float(
          original[aut].mean()-original[~aut].mean()),
      "L_island_bootstrap_95percentile_ci_pp":[float(v) for v in np.quantile(valid[:,0],[.025,.975])],
      "builtup_autonomous_minus_absent_pp":float(
           100*first[aut,5].mean()-100*first[~aut,5].mean()),
      "builtup_island_bootstrap_95percentile_ci_pp":[float(v) for v in np.quantile(valid[:,1],[.025,.975])],
      "ALAN_and_builtup_have_different_predictor_units_no_direct_numeric_comparison":True,
      "conditioned_on_unmeasured_biological_survey_detection_and_urban_habitat":True,
      "same_source_islands_resampled_jointly_across_species":True
    }

def run(sites_file,globi,groups,ledger,annual_worldpop,climate_root,ghsl_source,out):
    require(not out.exists(),"Refuse overwrite source-verified original analysis")
    original=load_sites(sites_file)
    flora=original_accepted_species(groups,original,load_labels(globi))
    trait=direct_plant_traits(ledger)
    taxa=pd.DataFrame({"species":flora.species.drop_duplicates()})
    taxa["autonomy"]=taxa.species.map(trait.autonomous_selfing_capacity)
    _,_,_,selected,s0=prepare_support(original,flora,taxa,span=.02)
    require(len(selected)==703 and s0["strict_selfing_species"]==703,
            "Frozen 703 independently direct-autonomous species source drift")
    labels=taxa.set_index("species").autonomy.loc[sorted(selected)]
    annual,ann_meta=source_annual_population(original,annual_worldpop)
    require(ann_meta["source_population_valid_sites"]==69525,
            "Previously verified full 5year WorldPop source changed")
    eligible,pmeta=within_year_island_photo_design(annual)
    require(len(eligible)==64966 and eligible.island.nunique()==202,
            "Fixed full-WorldPop population cohort changed")
    climate,cmeta=matching_climate_gradient(eligible,climate_root)
    require(len(climate)==64944,"Fixed original WorldClim source expected 64,944 sites")
    built,bmeta=authentic_ghsl_source(ghsl_source,climate)
    pos,iy,sup=candidate_species_photo_support(built,flora,labels)
    fits=[];blocks=[];failed={}
    for j,sp in enumerate(sorted(selected)):
        if sp not in pos or sp not in iy:
            reason="NO_SOURCE_DIEL_INFORMATIVE_PHOTO_STRATA"
            failed[reason]=failed.get(reason,0)+1
            continue
        ans,reason=fit_source_species(built,pos[sp],iy[sp])
        if reason:
            failed[reason]=failed.get(reason,0)+1
            continue
        q,islands,A,B=ans
        aut=str(labels[sp])
        fits.append({**q,"autonomy":aut})
        blocks.append({"q":q,"autonomy":aut,
                       "islands":islands,"A":A,"B":B})
        if (j+1)%100==0:
            print("ALAN_GHSL_ORIGINAL_FIT_PROGRESS",j+1,flush=True)
    require(len(fits)>0,"No GHSL light+population climate original model fit")
    summary=groups_summary(fits)
    b2015=island_bootstrap(blocks,built.island.unique(),epoch="ghsl2015")
    b2020=island_bootstrap(blocks,built.island.unique(),epoch="ghsl2020")
    data={
       "schema":"ALAN_source_703_plant_species_prior2015_GHSL_R2023A_worldpop_worldclim_20261116_v1",
       "status":"INDEPENDENT_SOURCE_JRC_GHSL_URBAN_BUILTUP_COVARIATE_EXPLORATORY_NOT_CAUSAL",
       "original_source":{
         "source_plant_photo_site_years":77647,
         "source_islands":562,
         "source_plant_species_site_years":2012925,
         "source_703_independent_direct_autonomous_selfing_species":703,
         "source_WorldPop_sameyear_direct_ALAN_2019_2025":69525,
         "source_original_climate_2015_2020_not_temporal_annual_weather":True},
       "worldpop_complete_source":ann_meta,
       "photo_light_population_source":pmeta,
       "worldclim_photo_source":{
          "joint_eligible_photo_site_years":len(climate),
          "independent_islands":int(climate.island.nunique()),
          "climate_1970_2000_static_normal":True},
       "builtup_source_metadata":bmeta,
       "species_support":{
          **sup,
          "fixed_original_703_species":703,
          "source_GHSL_builtup_species_model_fitted":len(fits),
          "autonomous_fitted":int(sum(x["autonomy"]=="autonomous" for x in fits)),
          "absent_fitted":int(sum(x["autonomy"]=="absent" for x in fits)),
          "excluded_source_reasons":failed},
       "same_source_model_slope_differences":summary,
       "island_bootstrap_prior2015":b2015,
       "island_bootstrap_sensitivity2020":b2020,
       "limits":{
         "2015_builtup_prior_all_photo_observations":True,
         "2020_builtup_post_2019_photographs_sensitivity_only":True,
         "GHS_BUILT_S_estimated_satellite_built_area_not_ecological_habitat":True,
         "population_collinearity_with_light_unresolved":True,
         "photo_nonrecord_is_not_true_absence":True,
         "no_local_population_decline_or_fruit_seed_observed":True,
         "no_species_efficiency_or_actual_moth_pollen_delivery_tested":True,
         "human_lighting_source_not_measured_with_local_lux":True,
         "phylogenetic_residual_and_observer_bias_not_modeled":True,
         "no_restricted_species_location_source_rows_exported":True,
         "no_causal_ALAN_claim":True
       }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_GHSL_BUILTUP_703_SPECIES",json.dumps({
       "build":bmeta,"species":data["species_support"],
       "traits":summary,"boot2015":b2015,"boot2020":b2020
    },ensure_ascii=False,allow_nan=False),flush=True)
    return data

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("sites","globi","groups","ledger","annual_worldpop",
                "climate_root","ghsl_source","out"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.annual_worldpop,
        a.climate_root,a.ghsl_source,a.out)
