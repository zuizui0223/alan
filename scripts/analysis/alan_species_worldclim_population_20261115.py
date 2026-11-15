#!/usr/bin/env python3
"""ALAN 703 measured-autonomy plants: source-matched WorldPop plus WorldClim.

Fixed original photographic response: source-accepted plant species recorded
versus NOT photographed in another photographed botanical source cell.
Nonrecord is NOT independently surveyed plant absence, reproduction or decline.

For each original species, fit on identical authentic source photo-cell-years:
  A = prior-PCNL + photo effort + same-year exact-pixel population
  B = A + WorldClim 2.1 static BIO1 temperature + log1p(BIO12 annual rain)
with island×observation-year FE and equal-island weighting. BIO15 seasonality
is audited but not added post-hoc (two climate axes fixed before source result).

This is NOT a causal ALAN model; population and climate do not adjust habitat,
plant observation selection or at-flower night illuminance.
"""
from __future__ import annotations
import argparse, json, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import t as student_t
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY,load_sites,load_labels,original_accepted_species,require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import (
    prepare_support,percentile
)
from scripts.analysis.alan_mois_worldpop_exact_pixel_support_20261115 import join
from scripts.analysis.alan_species_worldpop_conditional_20261115 import (
    within_year_island_photo_design,candidate_species_photo_support,fit_species
)
from scripts.analysis.alan_worldclim_original_climate_20261115 import (
    sample_sites
)

N_BOOT=400
SEED=20261115
BASE=("loglight","effort","logpop")
CLIMATE=("bio1","logprec")
ALL=BASE+CLIMATE

def matching_climate_gradient(source,climate_folder):
    d,rmeta=sample_sites(source,climate_folder)
    d=d.loc[d.climate_complete].copy()
    d["logprec"]=np.log1p(d.bio12.astype(float))
    # No ecological or outcome-driven threshold: reapply only exposure and
    # survey opportunity eligibility after authentic climate NODATA removal.
    group=d.groupby(["island","year"],observed=True)
    g=group.loglight.agg(["size","min","max","median"])
    good=g.loc[g["size"].ge(6)&(g["max"]-g["min"]).ge(.02)]
    d=d.merge(good[["median"]].reset_index().rename(columns={"median":"cmid"}),
              on=["island","year"],how="inner",validate="many_to_one")
    d["chigh"]=d.loglight.gt(d.cmid).astype(int)
    both=d.groupby(["island","year","chigh"],observed=True).size().unstack(
        "chigh",fill_value=0)
    for i in (0,1):
        if i not in both:both[i]=0
    iy=both.loc[both[0].ge(2)&both[1].ge(2)].reset_index()[["island","year"]]
    d=d.merge(iy,on=["island","year"],validate="many_to_one")
    d=d.reset_index(drop=True)
    d["sid"]=np.arange(len(d),dtype=np.int32)
    d["iy"],_=pd.factorize(d.island.astype(str).str.cat(
        d.year.astype(str),sep="|"),sort=True)
    require(len(d)>0 and d.island.nunique()>10,
            "WorldClim matched population light contrast missing")
    z=d[list(ALL)]-d.groupby("iy",observed=True)[list(ALL)].transform("mean")
    zstd=z.std()
    return d,{
        "original_population_light_eligible_site_years":len(source),
        "source_climate_complete_site_years":int(source.shape[0]-(~source.climate_complete).sum())
                  if "climate_complete" in source else None,
        "population_worldclim_photo_site_years_after_new_light_gate":len(d),
        "population_worldclim_supported_islands":int(d.island.nunique()),
        "population_worldclim_island_years":int(d.iy.nunique()),
        "within_island_year_global_residual_predictor_corr":(
             {a:{b:float(z[a].corr(z[b])) for b in ALL} for a in ALL}),
        "within_island_year_predictor_sd":{k:float(v) for k,v in zstd.items()},
        "raster_origin_verified_per_BIO":rmeta,
        "climate_is_1970_2000_static_10arcminute_normal":True,
        "no_missing_climate_imputation":True,
    }

def fit_five(d,positive,informative_iy):
    s=d.loc[d.iy.isin(informative_iy)].copy()
    if len(s)<30:return None,"TOO_FEW_SITES"
    y=np.isin(s.sid.to_numpy(np.int64),positive,assume_unique=True).astype(float)
    G=s.island.nunique();K=s.iy.nunique();Np=int(y.sum())
    if G<2 or K<4 or Np<12:
        return None,"INSUFFICIENT_INFORMATIVE_SOURCE_SUPPORT"
    if Np>=len(s):
        return None,"NONRECORD_OPPORTUNITY_MISSING"
    cols=s[list(ALL)].copy()
    cols["outcome"]=y
    centered=(cols-cols.groupby(s.iy).transform("mean")).to_numpy(float)
    X=centered[:,:5]
    yy=centered[:,-1]
    island_ids,unique=pd.factorize(s.island,sort=True)
    w=1/np.bincount(island_ids)[island_ids].astype(float)
    scales=np.sqrt(np.sum(w[:,None]*X**2,axis=0)/w.sum())
    if not np.isfinite(scales).all() or np.min(scales)<=1e-8:
        return None,"NO_WITHIN_ISLAND_YEAR_CLIMATE_VARIATION"
    Z=X/scales
    fits=[]
    for p in (3,5):
        XX=Z[:,:p]
        xtwx=(XX.T*w)@XX
        if np.linalg.matrix_rank(xtwx)<p or np.linalg.cond(xtwx)>1e8:
            return None,"ORIGINAL_CLIMATE_LIGHT_POP_RANK_OR_COLLINEAR"
        bb=np.linalg.solve(xtwx,XX.T@(w*yy))
        raw=bb/scales[:p]
        if not np.isfinite(raw).all():
            return None,"NONFINITE_SOURCE_COEFFICIENT"
        fits.append(raw)
    base,full=fits
    # Exactly identical species survey/photo/population frame, only climate
    # covariates differ. Independent rank-check against previous A source code.
    orig,reason=fit_species(d,positive,informative_iy)
    if reason is not None:return None,"BASE_MODEL_NOT_REPRODUCED"
    o,_,_,_=orig
    require(abs(100*base[0]-o["beta_L_conditional_population_pp"])<1.e-7,
            "WorldClim model baseline does not reproduce original matched-pop source")
    q={
       "n_islands":int(G),"n_strata":int(K),
       "n_positive_source_photo_cells":Np,"n_source_site_years":len(s),
       "beta_L_no_climate_pp":float(100*base[0]),
       "beta_L_with_climate_pp":float(100*full[0]),
       "beta_pop_no_climate_pp":float(100*base[2]),
       "beta_pop_with_climate_pp":float(100*full[2]),
       "beta_bio1_conditional_pp_per_C":float(100*full[3]),
       "beta_logprec_conditional_pp":float(100*full[4]),
       "beta_L_change_adding_climate_pp":float(100*(full[0]-base[0]))
    }
    # Only island CR1 uncertainty for adequate independent island support.
    if G>=10:
        residual=yy-Z@(full*scales)
        inv=np.linalg.inv((Z.T*w)@Z)
        scores=np.zeros((G,5))
        np.add.at(scores,island_ids,Z*(w*residual)[:,None])
        dof=len(s)-K-5
        if dof>0:
            cov=((G/(G-1))*((len(s)-1)/dof))*inv@(scores.T@scores)@inv
            se=np.sqrt(np.maximum(np.diag(cov),0))/scales
            crit=student_t.ppf(.975,G-1)
            q["light_se_climate_pp"]=float(100*se[0])
            q["light_ci95_climate_pp"]=[
                 float(100*(full[0]-crit*se[0])),float(100*(full[0]+crit*se[0]))]
    # Sufficient statistics for TRUE whole-island joint first-stage re-fits.
    aa=np.zeros((G,5,5),float)
    bb=np.zeros((G,5),float)
    for i in range(G):
        ids=np.flatnonzero(island_ids==i)
        part=X[ids]
        aa[i]=(part.T*w[ids])@part
        bb[i]=part.T@(w[ids]*yy[ids])
    return (q,unique.astype(str).to_numpy(),aa,bb),None

def output_summary(rows):
    d=pd.DataFrame([{k:v for k,v in x.items() if k!="autonomy"} for x in rows])
    a=pd.Series([x["autonomy"] for x in rows])
    results={"n_fitted":len(d)}
    fields=("beta_L_no_climate_pp","beta_L_with_climate_pp",
            "beta_pop_no_climate_pp","beta_pop_with_climate_pp",
            "beta_L_change_adding_climate_pp")
    for name,mask in (("all",np.ones(len(a),dtype=bool)),
                      ("autonomous",a.eq("autonomous").to_numpy()),
                      ("absent",a.eq("absent").to_numpy())):
        z=d.loc[mask]
        results[name]={"species":len(z),
            "species_ge10_islands":int(z.n_islands.ge(10).sum()),
            "species_ge30_islands":int(z.n_islands.ge(30).sum()),
            **{f:percentile(z[f]) for f in fields}}
    for field in fields:
        yes=d.loc[a.eq("autonomous").to_numpy(),field]
        no=d.loc[a.eq("absent").to_numpy(),field]
        results[field+"_autonomous_minus_absent_pp"]=(
            float(yes.mean()-no.mean()) if len(yes) and len(no) else None)
    return results

def boot_islands(models,universe):
    kept=[x for x in models if x["q"]["n_islands"]>=10]
    auto=sum(x["autonomy"]=="autonomous" for x in kept)
    absent=len(kept)-auto
    if auto<25 or absent<25:
        return {"status":"INSUFFICIENT_INDEPENDENT_PLANT_CLASS_SUPPORT",
                "species":len(kept),"autonomous":auto,"absent":absent}
    n=len(kept);islands=list(universe);m=len(islands)
    lookup={x:i for i,x in enumerate(islands)}
    A=np.zeros((n,m,5,5))
    B=np.zeros((n,m,5))
    group=np.zeros(n,bool)
    for i,x in enumerate(kept):
        group[i]=x["autonomy"]=="autonomous"
        ix=np.array([lookup[k] for k in x["islands"]],dtype=int)
        A[i,ix,:,:]=x["AA"]
        B[i,ix,:]=x["BB"]
    point=np.linalg.solve(A.sum(axis=1),B.sum(axis=1)[...,None])[...,0]
    expected=np.array([x["q"]["beta_L_with_climate_pp"] for x in kept])
    require(abs(100*point[:,0]-expected).max()<1e-7,
            "Whole-island source climate stats cannot reproduce original fit")
    rng=np.random.default_rng(20261115)
    bootstrap=np.full(N_BOOT,np.nan)
    valid=0
    for i in range(N_BOOT):
        mul=np.bincount(rng.integers(0,m,m),minlength=m).astype(float)
        aa=np.einsum("g,ngjk->njk",mul,A,optimize=True)
        bb=np.einsum("g,ngj->nj",mul,B,optimize=True)
        try:
            b=np.linalg.solve(aa,bb[...,None])[...,0]
        except np.linalg.LinAlgError:continue
        if not np.isfinite(b).all():continue
        bootstrap[i]=100*(b[group,0].mean()-b[~group,0].mean())
    x=bootstrap[np.isfinite(bootstrap)]
    if len(x)<int(.8*N_BOOT):
        return {"status":"NOT_ESTIMABLE_SOURCE_CLIMATE_BOOTSTRAP_RANK",
                "n_valid_draws":len(x),"n_requested":N_BOOT}
    return {
      "status":"ORIGINAL_SOURCE_WHOLE_ISLAND_BOOTSTRAP_CLIMATE_PLUS_POPULATION",
      "ge10_species":n,"autonomous":auto,"absent":absent,
      "n_matched_islands":m,
      "n_bootstrap_requested":N_BOOT,
      "n_bootstrap_valid":len(x),
      "conditional_light_autonomous_minus_absent_pp":float(
              expected[group].mean()-expected[~group].mean()),
      "bootstrap_percentile_95_ci_pp":[float(a) for a in np.quantile(x,[.025,.975])],
      "island_resampling_conditional_on_MOIS_photo_pop_coverage":True,
    }

def run(site_file,globi,groups,ledger,pop_sources,climate_root,out):
    require(not out.exists(),"Source output overwrite refused")
    original=load_sites(site_file)
    flora=original_accepted_species(groups,original,load_labels(globi))
    trait=direct_plant_traits(ledger)
    plants=pd.DataFrame({"species":flora.species.drop_duplicates()})
    plants["autonomy"]=plants.species.map(trait.autonomous_selfing_capacity)
    _,_,_,fixed703,s0=prepare_support(original,flora,plants,span=.02)
    require((s0["strict_selfing_species"],len(fixed703))==(703,703),
            "Locked original source 703 changed")
    original_labels=plants.set_index("species").autonomy.loc[sorted(fixed703)]
    join_receipt,pop_all=join(original,pop_sources,return_data=True)
    require(join_receipt["matched_usable_site_years"]==23018,
            "Signed original same-pixel-population support changed")
    base,baseinfo=within_year_island_photo_design(pop_all)
    require(len(base)==20582 and base.island.nunique()==100,
            "Source-matched prior population model reference drift")
    p,climinfo=matching_climate_gradient(base,climate_root)
    site_ids,strata,cover=candidate_species_photo_support(
        p,flora,original_labels)
    fits=[];blocks=[];fail={}
    for j,sp in enumerate(sorted(fixed703)):
        if sp not in site_ids or sp not in strata:
            fail["NO_INFORMATIVE_CLIMATE_SITE_YEAR"]=fail.get(
                "NO_INFORMATIVE_CLIMATE_SITE_YEAR",0)+1
            continue
        ans,why=fit_five(p,site_ids[sp],strata[sp])
        if why is not None:
            fail[why]=fail.get(why,0)+1
            continue
        fit,islands,aa,bb=ans
        typ=str(original_labels[sp])
        fits.append({**fit,"autonomy":typ})
        blocks.append({"q":fit,"autonomy":typ,
                       "islands":islands,"AA":aa,"BB":bb})
        if (j+1)%100==0:
            print("ALAN_AUTHENTIC_POP_WORLDCLIM_SPECIES_FIT",j+1,flush=True)
    require(len(fits)>0,"No source-supported combined population/climate fit")
    outdata={
       "schema":"ALAN_703_plant_original_photoPCNL_sameyearWorldPop_WorldClim10m_static_v1",
       "status":"ORIGINAL_SOURCE_CLIMATE_POPULATION_CONDITIONAL_MODEL_EXPLORATORY_NOT_CAUSAL",
       "provenance":{
           "original_source_77k_photo_years":77647,
           "original_703_direct_measured_autonomy_species":703,
           "authenticated_pop_23k_exactyear_1km_pixel":True,
           "original_10m_bio_archive_SHA256":"00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf",
           "WorldClim_normal_period_1970_2000":True
       },
       "population_join":{k:join_receipt[k] for k in
          ("matched_usable_site_years","matched_distinct_islands",
           "source_zero_population_site_years","same_worldpop_pixel_required")},
       "original_population_light_cohort":baseinfo,
       "climate_matched_source":climinfo,
       "species_support":{**cover,"fixed_original_703":703,
                          "n_pop_climate_model_fitted":len(fits),
                          "fit_autonomous":sum(x["autonomy"]=="autonomous" for x in fits),
                          "fit_absent":sum(x["autonomy"]=="absent" for x in fits),
                          "nonfit_reasons":fail},
       "same_cohort_model_comparison":output_summary(fits),
       "bootstrap":boot_islands(blocks,p.island.unique()),
       "limitations":{
           "no_true_absence_or_plant_reproduction":True,
           "pop_photo_matched_only_nonrandom_30_percent_of_original":True,
           "worldclim_10m_static_normal_coarse_resolution":True,
           "no_observed_2019_2025_annual_or_monthly_actual_climate":True,
           "no_land_use_or_urbanization_causal_independence":True,
           "strong_PCNL_pop_collinearity_remains":True,
           "no_phylogeny_or_continental_occurrence":True,
           "not_reconfirming_prior_positive_photo_colour_or_pollen_limit_results":True,
           "no_raw_source_botanical_trait_by_location_exported":True
       }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(outdata,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_WORLDPOP_WORLDCLIM_JOINT_SOURCE_REAL",json.dumps({
      "climate_design":climinfo,
      "species_support":outdata["species_support"],
      "same_cohort_model_comparison":outdata["same_cohort_model_comparison"],
      "bootstrap":outdata["bootstrap"]
    },ensure_ascii=False,allow_nan=False),flush=True)
    return outdata

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("sites","globi","groups","ledger","pop_sources","climate_root","out"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.pop_sources,a.climate_root,a.out)
