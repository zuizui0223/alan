#!/usr/bin/env python3
"""ALAN tier08 actual precolour photo opportunity × published moth-evidence tiers.

Photo samples are genuine original-photo IDs from the *already frozen*
24,576 tier08 source of the 93.2M photo index. Match species by exact WCVP
accepted plant name from signed source (not lookalike floral morphology).
Independent published GloBI source states document flower visitation or a
database pollinates claim, NOT a definitive local nocturnal, effective or
exclusive moth pollination. Walton 2020 documented nocturnal pollen carriage
without proving local stigma deposition, seed set or moth dependence.

Do NOT alter any of the 48 concurrently running original RGB worker photo
rosters because outcome-dependent sampling would destroy held-out comparability.
Report primary moth-associated subgroup, moth-database pollination-positive
subset, Walton field nocturnal pollen-transport evidence subset, and positive
butterfly comparison flag (not negative moth control). Source-only evidence
does not claim actual photographed RGB success, island-light effects, or power.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="ALAN_original24576_island_GloBI_moth_pollination_evidence_source_support_v1"
SHA={
 "queue":"5578cd2e02330ee951539d4b7d0a2e166b1c8b1ec6dc87f383f8625c7608098d",
 "accepted":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
 "globi":"a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994",
}
FLAGS=("moth_record_positive","moth_database_pollinates",
       "walton_direct_nocturnal_pollen_carriage","butterfly_record_positive")

def require(ok,reason):
    if not bool(ok):raise ValueError("ALAN_MOTH_REAL_SOURCE_EVIDENCE: "+reason)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def evidence_flags(globi):
    needed={"plant_name","guild","source_family","relation","study_title"}
    require(needed<=set(globi.columns),"Independent floral-interaction source fields missing")
    assert set(globi.guild).issubset({"moth_family_visit","butterfly_family_visit"})
    g=globi.copy()
    m=g.guild.eq("moth_family_visit")
    g["moth_record_positive"]=m
    g["moth_database_pollinates"]=m & g.relation.eq("pollinates")
    # The exact DOI/primary source link; a secondary DoPI study title citing
    # Walton is not substituted for original measured nocturnal collection.
    g["walton_direct_nocturnal_pollen_carriage"]=(
        m & g.study_title.str.contains("rsbl.2019.0877",regex=False,case=False,na=False)
    )
    g["butterfly_record_positive"]=g.guild.eq("butterfly_family_visit")
    flags=g.groupby("plant_name",observed=True)[list(FLAGS)].any().reset_index()
    require(flags.plant_name.is_unique,
            "No unique positive evidence label per published plant binomial")
    return flags

def join_original_queue(queue, accepted, globi):
    require(queue.photo_id.is_unique and queue.observation_id.is_unique and
            len(queue)==24576 and queue.wcvp_accepted_plant_name_id.ne("").all(),
            "24,576 source-original photo UUIDs or accepted taxa missing")
    require(queue.wave.value_counts().to_dict()=={"1":8192,"2":8192,"3":8192},
            "Frozen three 8192 original-photo source waves changed")
    names=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(names.wcvp_accepted_plant_name_id.is_unique,
            "Official WCVP accepted name ID maps to distinct botanical names")
    q=queue.merge(names,on="wcvp_accepted_plant_name_id",
                  how="left",validate="many_to_one")
    require(q.wcvp_accepted_name.notna().all() and
            len(q)==len(queue),
            "Some selected original photos lack WCVP accepted names")
    flags=evidence_flags(globi)
    q=q.merge(flags,left_on="wcvp_accepted_name",right_on="plant_name",
              how="left",validate="many_to_one")
    require(len(q)==len(queue) and q.photo_id.is_unique,
            "Published plant-guild source multiplied the original photo rows")
    for key in FLAGS:
        q[key]=q[key].fillna(False).astype(bool)
    return q,flags

def support_for(data,name):
    d=data.loc[data[name]].copy()
    site=d.groupby("target_id",observed=True).agg(
        n=("photo_id","size"),species=("wcvp_accepted_plant_name_id","nunique"),
        island=("land_id","first"))
    potential=site.loc[site.n.ge(5)&site.species.ge(3)]
    by_island=d.groupby("land_id",observed=True).target_id.nunique()
    mixed=d.assign(light=pd.to_numeric(d.log_light,errors="coerce")).groupby(
         "land_id",observed=True).light.agg(["min","max"])
    adequate=(mixed["max"]-mixed["min"]).ge(.05)
    return {
       "original_photo_attempts_not_ROI_success":int(len(d)),
       "recorded_plant_species":int(d.wcvp_accepted_plant_name_id.nunique()),
       "distinct_islands_with_any_photo":int(d.land_id.nunique()),
       "distinct_source_site_years":int(d.target_id.nunique()),
       "islands_with_at_least_2_source_photo_sites":int(by_island.ge(2).sum()),
       "islands_with_at_least_2_observed_sites_and_0p05_logVIIRS_gap":
           int((by_island.ge(2)&adequate.reindex(by_island.index,fill_value=False)).sum()),
       "precolour_site_years_eligible_for_5_attempt_3_species_minimum":
           int(len(potential)),
       "precolour_islands_with_any_such_5attempt_3species_site":
           int(potential.island.nunique()),
       "source_original_photo_attempts_by_wave":{
           str(k):int(v) for k,v in d.wave.value_counts().items()},
       "source_original_photo_attempts_by_area_third":{
           str(k):int(v) for k,v in d.area_third.value_counts().items()},
       "not_true_nocturnal_effectiveness":True,
       "actual_photo_numeric_ROI_success_unknown_until_original_RGB_collector":True,
    }

def run(queue,accepted,globi,out):
    paths={"queue":Path(queue),"accepted":Path(accepted),"globi":Path(globi)}
    for k,path in paths.items():
        require(sha(path)==SHA[k],"SHA-pinned published original source bytes drift: "+k)
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of source-validated moth selection audit")
    q=pd.read_csv(queue,compression="gzip",dtype=str,keep_default_na=False)
    a=pd.read_csv(accepted,compression="gzip",dtype=str,
        usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    g=pd.read_csv(globi,compression="gzip",dtype=str)
    require(len(g)==61967 and len(q)==24576,"Signed original GloBI or tier08 source size changed")
    data,flags=join_original_queue(q,a,g)
    result={
      "schema":SCHEMA,
      "source_SHA256":SHA,
      "source_original_RGB_photo_attempts":int(len(data)),
      "source_GSHHG_islands":int(data.land_id.nunique()),
      "source_photo_siteyears":int(data.target_id.nunique()),
      "source_GloBI_signed_positive_species_by_tier":{
          flag:int(flags[flag].sum()) for flag in FLAGS},
      "tier_source_support":{
          flag:support_for(data,flag) for flag in FLAGS},
      "source_nightly_pollen_transport_is_not_local_effective_pollination":True,
      "recorded_moth_visitation_not_moth_exclusivity":True,
      "no_record_not_negative_moth_or_negative_butterfly_evidence":True,
      "three_predeclared_full_photo_RGB_waves_unaltered":True,
      "no_original_ROI9_inference_or_new_RGB_pixels":True,
      "no_verified_white_flower_native_flora_or_ALAN_causal_effect":True,
    }
    positive_without_butterfly=data.moth_record_positive&~data.butterfly_record_positive
    result["moth_source_positive_without_known_butterfly_record_photo_attempts"]=int(
          positive_without_butterfly.sum())
    result["moth_source_positive_without_butterfly_record_NOT_proven_moth_only"]=True
    out.mkdir(parents=True)
    (out/"alan_tier08_24576_moth_positive_original_photo_support.json").write_text(
        json.dumps(result,sort_keys=True,indent=2)+"\n")
    pd.DataFrame([{
        "moth_evidence_tier":label,**{
           key:val for key,val in evidence.items() if not isinstance(val,dict)}
        } for label,evidence in result["tier_source_support"].items()]).to_csv(
        out/"moth_evidence_photo_and_island_source_support.csv",index=False)
    print("ALAN_TIER08_24576_MOTH_SOURCE_INDEPENDENT_METADATA_REAL_EVIDENCE",
          json.dumps({k:result[k] for k in ("source_original_RGB_photo_attempts",
            "source_GSHHG_islands","source_photo_siteyears",
            "tier_source_support")},sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("queue","accepted","globi","out"):
        p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args()
    run(args.queue,args.accepted,args.globi,args.out)
if __name__=="__main__":main()
