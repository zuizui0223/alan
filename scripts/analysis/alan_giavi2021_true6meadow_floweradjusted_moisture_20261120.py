#!/usr/bin/env python3
"""Six-meadow original Giavi 2021 ALAN wetness visitor analysis with floral cover.

Original Zenodo 4540407: 1884 source records, 6 lit-dark meadow pairs,
37 paired visits, 21 plant labels, 3 insect guilds, observed source zeros.
Only source-observed C/T matched plant-meadow-visit records, no assigned
biological zeros. Source Kew F_PLANTATT and F_Doring are independent traits.
Outcome Δlog1p(sum daytime visitors across 3 guilds) by same plant and
same paired visit; covariate Δlog1p(author-measured flowering cover).
Pair fixed effects (37), 2^6 whole-meadow flips (descriptive sharp-null).
Six meadows, NOT 193 independent ALAN treatments. A wetness effect
requires >=3 meadow sites with >=2 wet and >=2 drier scored plants and
>=4 unique wet and >=4 unique drier plants; else source-support STOP.
"""
from __future__ import annotations
import argparse,collections,itertools,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_giavi2021_originalRdata_meadowpair_support_20261120 import fetch
from scripts.analysis.alan_giavi2021_plantatt_moisture_dayvisits_20261120 import source_name,add_doring_only_records
from scripts.analysis.alan_anic2022_13species_kew_plantatt_20261120 import pull_kew_binary,trait_lookup,RAW_GIT_BLOB
SCHEMA="alan_Giavi2021_original6_meadows_flower_adjusted_wetness_moderation_v1"
GUILDS=("Diptera","Hymenoptera","Coleoptera")
TRAITS=("F_PLANTATT","F_Doring")
def check(v,msg):
    if not v:raise ValueError("GIAVI6_ORIGINAL_SOURCE_STOP: "+msg)
def matched_source(obj):
    f=obj["DAT"].copy()
    check(len(f)==1884,"1884 original source rows")
    keys=["site","repeate","plant","treatment"]
    check(set(keys+["insect_group","dVis","plab","pair"]).issubset(f),"source schema drift")
    for k in keys+["insect_group","pair"]:f[k]=f[k].astype(str).str.strip()
    for k in ("dVis","plab"):f[k]=pd.to_numeric(f[k],errors="coerce")
    check(not f[["dVis","plab"]].isna().any().any(),"source unknown visits or flowers")
    check(f[["dVis","plab"]].ge(0).all().all(),"source negative observations")
    check(set(f.treatment)=={"C","T"} and set(f.insect_group)==set(GUILDS),"author experimental arms/guilds drift")
    check(f.site.nunique()==6 and len(f[["site","repeate"]].drop_duplicates())==37,"six meadow 37 sampled pair source")
    check(len(f[["site","repeate","pair"]].drop_duplicates())==37,"original pair identity ambiguity")
    gr=f.groupby(keys,sort=True)
    check(gr.size().eq(3).all() and gr.insect_group.nunique().eq(3).all(),
          "each observed source plant-arm must contain 3 guild entries")
    check(gr.plab.nunique().eq(1).all(),"source flowering counts differ among 3 guild entries")
    g=gr.agg(visits=("dVis","sum"),flowers=("plab","first")).reset_index()
    for guild in GUILDS:
        sub=f.loc[f.insect_group==guild,keys+["dVis"]].rename(columns={"dVis":"v_"+guild})
        check(not sub.duplicated(keys).any(),"original duplicate guild row")
        g=g.merge(sub,on=keys,how="left",validate="one_to_one")
    check(not g.isna().any().any(),"original species-guild source merge failure")
    wide=g.pivot(index=["site","repeate","plant"],columns="treatment",
                 values=["visits","flowers",*[("v_"+x) for x in GUILDS]])
    paired=wide.dropna(subset=[("visits","C"),("visits","T"),("flowers","C"),("flowers","T")])
    check(len(paired)==193,"source matched plant x site x sampling event drift")
    d=paired.reset_index()
    d.columns=[a if not b else (a+"_"+b) for a,b in d.columns]
    d["pair_key"]=d["site"]+"|"+d["repeate"]
    d["delta_flower"]=np.log1p(d["flowers_T"])-np.log1p(d["flowers_C"])
    for col in ("visits",*("v_"+x for x in GUILDS)):
        d["delta_"+col]=np.log1p(d[col+"_T"])-np.log1p(d[col+"_C"])
    return d,{"rows":len(f),"meadow_pairs":f.site.nunique(),
              "sampling_events":len(f[["site","repeate"]].drop_duplicates()),
              "source_plant_arm_cells":len(g),"matched_plant_pair_cells":len(d),
              "observed_source_zero_visits":int(f.dVis.eq(0).sum()),
              "source_guild_flower_cover_identical":True}
def add_moisture(d):
    raw=pull_kew_binary()
    t,n=trait_lookup(raw);t=add_doring_only_records(raw,t)
    out=d.copy();names={}
    for name in sorted(d.plant.unique()):
        scientific=source_name(name)
        source=t.get(scientific) if scientific else None
        names[name]={"binomial":scientific,
          "F_PLANTATT":float(source["F_PLANTATT"]) if source and source["F_PLANTATT"] is not None else None,
          "F_Doring":float(source["F_Doring"]) if source and source["F_Doring"] not in ("NA",None,"") else None}
    for k in TRAITS:out[k]=out.plant.map({z:x[k] for z,x in names.items()})
    return out,names,n
def source_support(z,trait):
    plants=z[["plant",trait]].drop_duplicates()
    wet=int(plants[trait].gt(7).sum());dry=int(plants[trait].le(7).sum())
    site={}
    for s,g in z.groupby("site"):
        pairs=g[["plant",trait]].drop_duplicates()
        site[s]={"wet":int(pairs[trait].gt(7).sum()),
                 "drier":int(pairs[trait].le(7).sum()),
                 "matched_plant_sampling_rows":len(g)}
    mixed=sum(x["wet"]>0 and x["drier"]>0 for x in site.values())
    rich=sum(x["wet"]>=2 and x["drier"]>=2 for x in site.values())
    return {"source_meadow_sites":len(site),"source_scored_species":len(plants),
        "source_wet_species_gt7":wet,"source_drier_species_le7":dry,
        "sites_with_wet_and_drier":mixed,"sites_with_2wet_2drier":rich,
        "site_counts":site,"confirmatory_moisture_moderator_support":wet>=4 and dry>=4 and rich>=3,
        "status":"SOURCE_WET_DRY_COMMON_SUPPORT_STOP" if not(wet>=4 and dry>=4 and rich>=3)
            else "COMMON_SUPPORT_PASS_REPLICATION_STILL_REQUIRED"}
def estimate(z,trait,response="visits",flower=True):
    if len(z)<8 or z.plant.nunique()<3:return None
    levels=sorted(z.pair_key.unique())
    fe=np.stack([z.pair_key.eq(k).to_numpy(float) for k in levels],axis=1)
    X=np.column_stack([fe,z[trait].to_numpy(float)-5.])
    if flower:X=np.column_stack([X,z.delta_flower.to_numpy(float)])
    y=z["delta_"+response].to_numpy(float)
    if not np.isfinite(X).all() or not np.isfinite(y).all() or np.linalg.matrix_rank(X)<X.shape[1]:return None
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    return {"beta_F_per_ordinal_point":float(beta[len(levels)]),
      "beta_delta_log1p_flower":float(beta[-1]) if flower else None,
      "source_matched_plant_pair_records":len(z),"source_true_site_pairs":z.site.nunique(),
      "unique_source_species":z.plant.nunique(),"sampling_pair_FE_count":len(levels),
      "model":"paired lit-minus-dark log1p visits with sampling-event fixed effects and optional measured floral cover",
      "no_causal_mediation_claim":True}
def site_signflip(z,trait,response="visits",flower=True):
    sites=sorted(z.site.unique())
    if len(sites)<2:return None
    original=estimate(z,trait,response,flower)
    if original is None:return None
    vals=[]
    for signs in itertools.product([-1.,1.],repeat=len(sites)):
        zz=z.copy();sgn=zz.site.map(dict(zip(sites,signs))).to_numpy(float)
        zz["delta_"+response]=sgn*z["delta_"+response].to_numpy(float)
        zz["delta_flower"]=sgn*z.delta_flower.to_numpy(float)
        res=estimate(zz,trait,response,flower)
        if res is not None:vals.append(res["beta_F_per_ordinal_point"])
    b=original["beta_F_per_ordinal_point"]
    return {"observed_beta":b,"exact_whole_meadow_label_flips":len(vals),
       "true_original_meadows":len(sites),
       "two_sided_abs_extreme_fraction":sum(abs(v)>=abs(b)-1e-10 for v in vals)/len(vals),
       "not_original_randomized_unit_test":True}
def run():
    dat,proof=matched_source(fetch())
    f,labels,n=add_moisture(dat)
    out={}
    for trait in TRAITS:
        full=f[f[trait].notna()].copy()
        pcount=full.groupby("plant").site.nunique()
        strict=full[full.plant.isin(pcount[pcount>=2].index)].copy()
        cross=[s for s,g in full.groupby("site") if
               g.loc[g[trait]>7,"plant"].nunique()>0 and
               g.loc[g[trait]<=7,"plant"].nunique()>0]
        mixed=full[full.site.isin(cross)].copy()
        cohorts={"all_scored_true_plant_pairs":full,"plants_observed_2plus_meadows":strict,
                 "same_meadow_wet_drier_overlap_only":mixed}
        co={}
        for key,z in cohorts.items():
            models={}
            for g in ("visits",*("v_"+x for x in GUILDS)):
                models[g]={"not_flower_adjusted":estimate(z,trait,g,False),
                           "flower_adjusted":estimate(z,trait,g,True)}
            # The plant floral-resource proxy is a separate observational
            # outcome, NOT a proven causal mediator of experimental ALAN.
            floral=estimate(z,trait,"flower",False)
            models["flower_cover"]={"lit_minus_dark_log1p_flower_as_response":floral}
            total_raw=models["visits"]["not_flower_adjusted"]
            total_adj=models["visits"]["flower_adjusted"]
            if floral and total_raw and total_adj:
                algebra=float(floral["beta_F_per_ordinal_point"]*
                              total_adj["beta_delta_log1p_flower"])
                differential=(total_raw["beta_F_per_ordinal_point"]-
                              total_adj["beta_F_per_ordinal_point"])
                check(abs(algebra-differential)<1e-8,
                      "original same-cohort OLS FWL floral decomposition not exact")
                models["flower_cover"]["FWL_association_partition"]={
                  "unadjusted_F_slope_minus_flower_adjusted_F_slope":differential,
                  "concurrent_flower_F_slope_times_visit_flower_association":algebra,
                  "within_source_fitted_on_identical_true_plant_pairs":True,
                  "this_is_algebra_not_causal_mediation":True}
            if key in ("all_scored_true_plant_pairs","same_meadow_wet_drier_overlap_only"):
                models["visits"]["site_signflip_flower_adjusted"]=site_signflip(z,trait)
                models["visits"]["site_signflip_not_flower_adjusted"]=site_signflip(z,trait,flower=False)
                models["flower_cover"]["exact_meadow_signflip_source_flower_F"]=(
                    site_signflip(z,trait,response="flower",flower=False))
                models["flower_cover"]["leave_one_true_meadow_out_F_flower"]={
                    site:(lambda r:r["beta_F_per_ordinal_point"] if r else None)(
                       estimate(z.loc[z.site!=site],trait,"flower",False))
                    for site in sorted(z.site.unique())}
                models["visits"]["leave_true_meadow_out_F_betas"]={
                    site:(lambda r:r["beta_F_per_ordinal_point"] if r else None)(
                        estimate(z.loc[z.site!=site],trait,"visits",True))
                    for site in sorted(z.site.unique())}
            co[key]={"paired_rows":len(z),"independent_original_meadow_pairs":z.site.nunique(),
                "named_species":z.plant.nunique(),"wet_species":z.loc[z[trait]>7,"plant"].nunique(),
                "drier_species":z.loc[z[trait]<=7,"plant"].nunique(),"guilds":models}
        out[trait]={"source_support":source_support(full,trait),"cohorts":co}
    return {"schema":SCHEMA,"original_Giavi2021_DOI":"10.1038/s41467-021-22011-8",
      "original_data_DOI":"10.5281/zenodo.4540407","original_Zenodo_RData_MD5":"f60b21d82c163cbb53520fda0bc61985",
      "Kew_trait_git_blob_SHA1":RAW_GIT_BLOB,"Kew_source_rows":n,
      "original_source":proof,"source_named_plant_traits":labels,"trait_results":out,
      "no_193_independent_randomized_treatments":True,
      "causal_mediation_from_concurrent_floral_cover_not_identified":True,
      "observed_daytime_visits_only_not_pollen_or_seed":True,
      "unobserved_plant_arm_not_assumed_zero":True,
      "cannot_establish_wet_plant_macro_mechanism":True}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);a=p.parse_args()
    v=run();target=Path(a.out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(v,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    report={}
    for t,d in v["trait_results"].items():
        report[t]={"support":{k:z for k,z in d["source_support"].items() if k!="site_counts"},
            "results":{q:{"matched_plant_pair_rows":z["paired_rows"],
                           "true_meadows":z["independent_original_meadow_pairs"],
                           "plants":z["named_species"],
                           "guild_models":z["guilds"]}
                       for q,z in d["cohorts"].items()}}
    print("GIAVI2021_ACTUAL_FLOWER_ADJUSTED_F_SOURCE_AND_6SITE_NULL "+
          json.dumps(report,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
