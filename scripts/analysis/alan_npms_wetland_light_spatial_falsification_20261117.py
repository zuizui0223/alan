#!/usr/bin/env python3
"""Falsification audit for authentic NPMS wetland × MOIST × prior PCNL positive-cover finding.

Scientific question: does the previously source-supported triple interaction
survive light time-series switches, independently documented brightening,
shared biological taxa, and coarse spatial-block dependence?

All source values are read from original published EIDC NPMS and the pinned
2015–25 source-positive plant-cover artifact. Individual plot/geocoordinates,
volunteer IDs, source biological rows NEVER leave RAM or appear in output.
No biological plant absences are imputed, no causal/mechanism mediation claims.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from scripts.acquisition import alan_npms_original_habitat_source_gate_20261117 as source
from scripts.analysis import alan_npms_independent_moisture_priorlight_20261117 as core
from scripts.analysis import alan_npms_true_wetland_moisture_priorlight_20261117 as habitat

SCHEMA="alan_npms_wetland_light_timeseries_spatial_falsification_v1"
GATE={"min_wetland_plots":30,"min_nonwetland_plots":100,
      "min_wetland_shared_years":75,"min_nonwetland_shared_years":200,
      "min_each_wetness_guild_species":15}
TERMS=["light_moist","light_moist_wetland","doy30"]
FE=["panel_key","species_year","plotyear_method_season"]

def require(cond, message):
    if not bool(cond):raise ValueError("ALAN_LIGHT_FALSIFICATION_STOP: "+message)

def common_support(d):
    wet=d.loc[d.inland_wetland]
    other=d.loc[~d.inland_wetland]
    def items(a):
        return {
            "plots":int(a.location_id.nunique()),
            "rows":int(len(a)),
            "distinct_taxa":int(a.species.nunique()),
            "moist_gt7_taxa":int(a.loc[a.wet7==1,"species"].nunique()),
            "moist_le7_taxa":int(a.loc[a.wet7==0,"species"].nunique()),
            "joint_source_plotyears":int(a.plotyear_method_season.nunique())
        }
    a,b=items(wet),items(other)
    valid=(a["plots"]>=GATE["min_wetland_plots"] and b["plots"]>=GATE["min_nonwetland_plots"]
        and a["joint_source_plotyears"]>=GATE["min_wetland_shared_years"]
        and b["joint_source_plotyears"]>=GATE["min_nonwetland_shared_years"]
        and min(a["moist_gt7_taxa"],a["moist_le7_taxa"],
                b["moist_gt7_taxa"],b["moist_le7_taxa"])>=GATE["min_each_wetness_guild_species"])
    return {"inland_wetland":a,"other":b,"pass":bool(valid)}

def mixed_reselect(data):
    # A plant-species subset may lose all high/low-MOIST co-observation in
    # a source site-year; re-gate rather than retain artificially shared cells.
    if data.empty:return data
    y=data.groupby("plotyear_method_season").agg(
        minimum=("wet7","min"),maximum=("wet7","max"),
        taxa=("species","nunique"))
    keep=y.index[(y.minimum==0)&(y.maximum==1)&(y.taxa>=2)]
    return data[data.plotyear_method_season.isin(keep)].copy()

def plot_light_timeseries(eligible,all_source):
    l=all_source[["location_id","year","prior_light"]].drop_duplicates()
    require(not l.duplicated(["location_id","year"]).any(),
            "multiple original satellite light records for same physical plot-year")
    wet_by_plot=eligible[["location_id","inland_wetland"]].drop_duplicates()
    require(not wet_by_plot.location_id.duplicated().any(),
            "a physical plot cannot be both wetland and another habitat")
    l=l.merge(wet_by_plot,on="location_id",how="inner",validate="many_to_one")
    site=[]
    for plot,g in l.groupby("location_id"):
        g=g.sort_values("year")
        light=g.prior_light.to_numpy(float)
        years=g.year.to_numpy(float)
        require(np.isfinite(light).all() and np.isfinite(years).all(),
                "missing original PCNL or year")
        n=len(years)
        valid_pair=n>=2
        slope=(float(np.polyfit(years-years.mean(),light,1)[0]) if valid_pair else None)
        adj=np.abs(np.diff(light)/np.diff(years)) if valid_pair else np.array([])
        site.append({
            "location_id":plot,
            "inland_wetland":bool(g.inland_wetland.iloc[0]),
            "source_observed_years":n,
            "light_has_zero_and_positive":bool((light<=1e-12).any() and (light>1e-12).any()),
            "source_light_range":float(np.ptp(light)),
            "large_annual_light_jump_1logunit":bool(np.any(adj>=1.)),
            "light_slope_per_year":slope,
        })
    out=pd.DataFrame(site)
    require(out.location_id.is_unique,"duplicate independent source plot")
    return out

def group_by_sign(s,mode):
    if mode=="brightening":
        return s[(s.light_slope_per_year>0)].location_id
    if mode=="darkening":
        return s[(s.light_slope_per_year<0)].location_id
    raise ValueError(mode)

def estimate(common,mode):
    support=common_support(common)
    result={"subset":mode,"support":support}
    if not support["pass"]:
        result["status"]="original_wetland_common_support_STOP"
        return result
    fit=core.fit(common,TERMS,FE,two_way=True)
    result["status"]=fit["status"]
    if fit["status"]!="observational_positive_cover_only":
        result["fit_status"]=fit["status"]
        return result
    e=fit["coefficients"]["light_moist_wetland"]
    result["triple_habitat_moisture_light"]={
        "beta":e["beta"],"plot_cluster_SE":e["SE_real_plot_cluster"],
        "plot_95CI":e["nominal_95CI"],"plot_nominal_p":e["nominal_two_sided_p"],
        "two_way":fit.get("two_way_plot_species_robustness")
    }
    return result

def source_spatial_blocks(npms,common):
    # Exact source co-ordinates are used transiently to assign large 1-degree
    # geography groups. Neither coordinates nor plot IDs are emitted.
    original=source.read(npms)["sample"]
    xy=original[["location_id","LATITUDE","LONGITUDE"]].copy()
    xy["lat"]=pd.to_numeric(xy.LATITUDE,errors="coerce")
    xy["lon"]=pd.to_numeric(xy.LONGITUDE,errors="coerce")
    xy=xy.dropna(subset=["lat","lon"])
    xy=xy[xy.lat.between(47,62)&xy.lon.between(-13,5)]
    centers=xy.groupby("location_id")[["lat","lon"]].median()
    require(centers.index.is_unique,"duplicate coordinate center")
    centers["grid1"]=centers.apply(
        lambda r:f"{math.floor(r['lat'])}_{math.floor(r['lon'])}",axis=1)
    centers["grid2"]=centers.apply(
        lambda r:f"{math.floor(r['lat']/2)*2}_{math.floor(r['lon']/2)*2}",axis=1)
    common=common.join(centers[["grid1","grid2"]],on="location_id")
    require(common.grid1.notna().all(),"unlocated original plant plot source")
    return common

def clustered_grid_inference(common,grid):
    n=len(common)
    w=1/common.groupby("location_id").location_id.transform("size").to_numpy(float)
    Z,notes=core.within_residualize(common,"cover_pp",TERMS,FE,w)
    Y,X=Z[:,0],Z[:,1:]
    k=X.shape[1]
    M=(X*w[:,None]).T@X
    require(np.linalg.matrix_rank(M)==k,"within-site model not estimable")
    B=np.linalg.inv(M)
    beta=B@(X*w[:,None]).T@Y
    residual=Y-X@beta
    score=X*(w*residual)[:,None]
    def cov(labels):
        ix=pd.factorize(labels)[0]
        g=int(ix.max()+1)
        sums=np.zeros((g,k))
        np.add.at(sums,ix,score)
        meat=(sums.T@sums)*g/(g-1)*n/(n-k)
        return meat,g
    gram,g=cov(common[grid])
    group_wet=common[[grid,"inland_wetland"]].drop_duplicates()
    wet_blocks=group_wet.loc[group_wet.inland_wetland,grid].nunique()
    result={"grid":grid,"total_spatial_blocks":g,
            "inland_wetland_blocks":int(wet_blocks)}
    if g<10 or wet_blocks<5:
        result["status"]="insufficient_coarse_spatial_blocks_STOP"
        return result
    V=B@gram@B
    if V[1,1]<=0:
        result["status"]="nonpositive_spatial_covariance_STOP"
        return result
    se=math.sqrt(V[1,1])
    q=t.ppf(.975,g-1)
    result.update(
        status="source_coarse_spatial_cluster_estimated",
        beta=float(beta[1]),SE=float(se),
        nominal_95CI=[float(beta[1]-q*se),float(beta[1]+q*se)],
        nominal_p=float(2*t.sf(abs(beta[1]/se),g-1)),
        numerical_FWL_source_beta_crosscheck=float(beta[1]),
    )
    # A sensitivity using BOTH spatial blocks and taxon dependence.
    mb,gb=cov(common[grid])
    ms,gs=cov(common.species)
    comb=common[grid].astype(str)+"#"+common.species.astype(str)
    mi,gi=cov(comb)
    vt=B@(mb+ms-mi)@B
    if vt[1,1]>0 and np.isfinite(vt[1,1]):
        sx=math.sqrt(vt[1,1])
        d=min(gb,gs)-1
        q=t.ppf(.975,d)
        result["two_way_grid_and_species"]={
            "spatial_blocks":gb,"plant_taxa":gs,
            "spatial_x_species_intersections":gi,
            "SE":float(sx),
            "nominal_95CI":[float(beta[1]-q*sx),float(beta[1]+q*sx)],
            "nominal_p":float(2*t.sf(abs(beta[1]/sx),d)),
        }
    else:
        result["two_way_grid_and_species"]={"status":"nonpositive_covariance_STOP"}
    return result

def execute(npms,cover,receipt,eiv):
    base,common,eligible=habitat.study(npms,cover,receipt,eiv,return_source_panels=True)
    require(base["source_gate"]["pass"],"original source population failed existing gate")
    assert common_support(common)["inland_wetland"]["plots"]==47
    site=plot_light_timeseries(eligible,common)
    all_active_sites=set(common.location_id)
    site=site[site.location_id.isin(all_active_sites)].copy()
    require(len(site)==326,"original 326-site wetland source frame changed")
    report={
        "schema":SCHEMA,
        "source_CRC_and_prior_pcNL_hash_locked":True,
        "original_coobserved_habitat_gate":common_support(common),
        "causal_ALAN_effect":False,
        "assessed_response":"only original positive-recorded NPMS DOMIN vegetation cover",
        "true_plant_local_absence":False,
        "survival_selection_not_resolved":True,
        "falsification_posthoc":True,
        "source_light_site_diagnostics":{},
        "robustness":{},
        "GATE":GATE,
    }
    for wet,grp in site.groupby("inland_wetland"):
        ss=grp.light_slope_per_year
        report["source_light_site_diagnostics"]["inland_wetland" if wet else "other"]={
            "physical_sites":len(grp),
            "source_zero_positive_switch_sites":int(grp.light_has_zero_and_positive.sum()),
            "source_large_jump_ge1_per_year_sites":int(grp.large_annual_light_jump_1logunit.sum()),
            "source_brightening_slope_sites":int((ss>0).sum()),
            "source_darkening_slope_sites":int((ss<0).sum()),
            "source_flat_slope_sites":int((ss==0).sum()),
            "median_site_source_light_slope_per_year":float(ss.median()),
            "median_original_source_light_range":float(grp.source_light_range.median()),
            "source_light_range_ge0p1_sites":int((grp.source_light_range>=.1).sum()),
            "source_light_range_ge0p3_sites":int((grp.source_light_range>=.3).sum()),
        }
    def run_subset(name,allowed,selection=None):
        subset=common.loc[common.location_id.isin(set(allowed))].copy()
        if selection is not None:subset=subset.loc[subset.species.isin(selection)].copy()
        subset=mixed_reselect(subset)
        report["robustness"][name]=estimate(subset,name)

    run_subset("exclude_same_site_zero_positive_switch",
        site.loc[~site.light_has_zero_and_positive,"location_id"])
    run_subset("exclude_source_jumps_ge1_per_year",
        site.loc[~site.large_annual_light_jump_1logunit,"location_id"])
    for direction in ("brightening","darkening"):
        run_subset("site_source_trend_"+direction,group_by_sign(site,direction))
    wet_species=set(common.loc[common.inland_wetland,"species"])
    other_species=set(common.loc[~common.inland_wetland,"species"])
    shared=wet_species & other_species
    report["shared_identical_species_across_actual_habitats"]=len(shared)
    run_subset("identical_species_documented_in_both_habitats",
               site.location_id,shared)
    # Restore original published site coordinates transiently for the
    # 1-degree/2-degree spatial sensitivity; no coordinates are exported.
    common=source_spatial_blocks(npms,common)
    # Stronger ecological control: do not assume the *same species* follows
    # one national annual trajectory across wetland and other habitats.
    # Source taxon × habitat × year FE isolates between-plot exposure effects
    # within each habitat and is a posthoc falsifier, never a primary result.
    common["species_habitat_year"]=(
         common.species.astype(str)+"|"+common.inland_wetland.astype(int).astype(str)
         +"|"+common.year.astype(str))
    report["robustness"]["species_by_habitat_by_year_FE"]={
        "status":"posthoc_extra_year_FE",
        "source_support":common_support(common),
        "fit":core.fit(common,TERMS,
           ["panel_key","species_habitat_year","plotyear_method_season"],
           two_way=True)
    }
    # Negative temporal-order check: a PCNL value from the observation year
    # AFTER the botanical survey is only available when the physical plot
    # is genuinely visited again two calendar years later. This future-year
    # predictor is NOT an imputed or forward-filled satellite measurement.
    # If the future-light association is similar, we cannot interpret the
    # old prior-light association as temporally causal.
    exposure=eligible[["location_id","year","prior_light"]].drop_duplicates()
    require(not exposure.duplicated(["location_id","year"]).any(),
            "source repeated plot-year PCNL has two values")
    future=exposure.rename(columns={"prior_light":"future_year_light"}).copy()
    future["year"]=future["year"]-2
    lead=common.merge(future,on=["location_id","year"],how="left",
                      validate="many_to_one")
    lead=lead[lead.future_year_light.notna()].copy()
    lead=mixed_reselect(lead)
    report["robustness"]["future_illumination_year_placebo"]={
        "status":"source_future_year_exposure_support",
        "support":common_support(lead),
        "source_logic":"observed at botanical year t+2 maps to PCNL observed year t+1; no future value is imputed"
    }
    if report["robustness"]["future_illumination_year_placebo"]["support"]["pass"]:
        lead["future_light_moist"]=lead.future_year_light*lead.moist_centered
        lead["future_light_moist_wetland"]=(
            lead.future_light_moist*lead.inland_wetland.astype(float))
        pmodel=core.fit(lead,
            ["future_light_moist","future_light_moist_wetland","doy30"],FE,two_way=True)
        report["robustness"]["future_illumination_year_placebo"]["fit"]=pmodel
        # Joint past/future source PCNL contrast may be unidentifiable owing
        # to strong serial correlation; allow the core rank gate to STOP.
        jmodel=core.fit(lead,
            ["light_moist","light_moist_wetland",
             "future_light_moist","future_light_moist_wetland","doy30"],
             FE,two_way=True)
        report["robustness"]["future_illumination_year_placebo"]["joint_old_future_fit"]=jmodel
    else:
        report["robustness"]["future_illumination_year_placebo"]["status"]="insufficient_future_exposure_common_support_STOP"
    for grid in ("grid1","grid2"):
        c=common[[grid,"location_id","inland_wetland"]].drop_duplicates()
        count=c.groupby(grid).agg(
            total_plots=("location_id","nunique"),
            wetland_plots=("inland_wetland","sum")).sort_values(
                    ["wetland_plots","total_plots"],ascending=False)
        report["source_light_site_diagnostics"][f"{grid}_spatial_support"]={
            "total_blocks":len(count),
            "blocks_with_wetland_plots":int((count.wetland_plots>0).sum()),
            "largest_block_wetland_plot_share":(
                float(count.wetland_plots.max()/47)),
            "largest_two_block_wetland_plot_share":(
                float(count.wetland_plots.iloc[:2].sum()/47)),
        }
        report["robustness"][f"{grid}_spatial_cluster"]=clustered_grid_inference(common,grid)
        for i,omit in enumerate(count.index[:2]):
            run_subset(f"{grid}_exclude_top_wetland_plot_block_{i+1}",
                       common.loc[common[grid]!=omit,"location_id"])
    return report

def main():
    p=argparse.ArgumentParser()
    for name in ("npms","cover","receipt","eiv","out"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    result=execute(a.npms,a.cover,a.receipt,a.eiv)
    target=Path(a.out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_NPMS_TRUE_WETLAND_EXPOSURE_SPATIAL_FALSIFICATION",json.dumps(
        {"light_site_diagnostics":result["source_light_site_diagnostics"],
         "shared_species":result["shared_identical_species_across_actual_habitats"],
         "robustness":{k:{
             "status":v["status"],
             "wet_plots":v.get("support",{}).get("inland_wetland",{}).get("plots"),
             "beta":v.get("triple_habitat_moisture_light",{}).get("beta",
                    v.get("beta")),
             "CI":v.get("triple_habitat_moisture_light",{}).get("plot_95CI",
                  v.get("nominal_95CI"))}
             for k,v in result["robustness"].items()}}
        ,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
