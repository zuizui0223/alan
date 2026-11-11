#!/usr/bin/env python3
"""Gate additional tier09 moth-documented photos against ALL 112428 real measured-photo IDs.

Source-only: exactly the source frame already fitted to ALAN, with genuine
numeric original RGB ROI9 status. A GloBI-claim-positive source photograph is
not necessarily moth-dependent, nighttime pollinated, or a biological white
flower. Pre-measurement support is not observed ROI9 or causal ALAN evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

SCHEMA="ALAN_tier09_vs_entire_measured112428_RGB_source_geography_v1"


def require(ok, msg):
    if not bool(ok):
        raise ValueError("ALAN_TIER09_FULL_GEOGRAPHY: "+msg)


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()


def boolean(series, name):
    d=series.astype(str).str.lower().map({"true":True,"false":False})
    require(d.notna().all(), "nonbinary status in "+name)
    return d.astype(bool)


def evaluate(photo, cells, candidate, island_year, allmodel):
    require(len(photo)==112428 and photo.photo_id.is_unique,
            "pooled model exact source photo IDs not fully distinct")
    require(boolean(photo.real_numeric_ROI9,"pooled ROI9").sum()==40914,
            "pooled original photographic numeric ROI9 drift")
    require(candidate.photo_id.is_unique and len(candidate)==78214,
            "tier09 original-photo candidate ID roster changed")
    duplicates=set(photo.photo_id)&set(candidate.photo_id)
    require(not duplicates,"new source repeats previously attempted exact photo ID")
    require(allmodel["schema"]=="ALAN_true112428_original_photo_RGB_GloBI_moth_guild_stratified_v1",
            "full model schema drift")
    base=allmodel["summary"]["all_original_photographed_flora"]
    require(base["multisite_islands"]==133 and base["multisite_site_years"]==3315 and
            abs(base["island_equal_area_ALAN"]["beta"]-0.6031352440580511)<1e-4,
            "pre-registered full model geography/term drift")
    p=photo.copy()
    parts=p.target_id.str.split("|",regex=False,expand=True)
    require(parts.shape[1]==3 and
            parts[0].eq(p.year.astype(str)).all() and
            parts[2].eq(p.land_id.astype(str)).all(),
            "real original measured photo target is not year|0.05-cell|GSHHG")
    p["cell_id"]=parts[1]
    p["real_numeric_ROI9"]=boolean(p.real_numeric_ROI9,"pooled ROI9")
    c=cells.copy()
    c["strict_photo_gate"]=boolean(c.strict_photo_gate,"original strict site gate")
    cp=c.target_id.str.split("|",regex=False,expand=True)
    require(cp.shape[1]==3 and cp[0].eq(c.year.astype(str)).all() and
            cp[2].eq(c.land_id.astype(str)).all() and c.target_id.is_unique,
            "signed pooled accepted-taxon cell geography drift")
    c["cell_id"]=cp[1]
    require(int(c.strict9_success.sum())==40914,
            "source cell true ROI9 successes do not sum to signed original photos")
    n=candidate.copy()
    n["viirs_eligible"]=boolean(n.viirs_eligible,"existing-cell same-year VIIRS")
    keep=island_year.loc[island_year.cohort.eq("GloBI_moth_record_positive")].copy()
    require(len(keep)==137 and keep.land_id.nunique()==56 and
            keep[["land_id","year"]].drop_duplicates().shape[0]==137,
            "verified tier09 independent light-support island-year pairs drift")
    keep_ids=set(keep.land_id)
    require(keep_ids.issubset(set(n.land_id)),"qualified islands lack signed source")
    q=(n.loc[n.viirs_eligible & n.land_id.isin(keep_ids)]
       .groupby(["land_id","year","cell_id"],observed=True)
       .agg(photo_IDs=("photo_id","nunique"),
            accepted_taxa=("species","nunique"),
            area_km2=("area_km2","first")).reset_index())
    q=q.loc[q.photo_IDs.ge(5)&q.accepted_taxa.ge(3)].copy()
    require(len(q)==2328 and q[["land_id","year"]].drop_duplicates().shape[0]>=56,
            "previous source-only qualified tier09 light-island cell support drift")
    # Use original GSHHG island and same-year exact 0.05-degree cell equality.
    keycols=["land_id","year","cell_id"]
    def keys(df):
        return set(map(tuple,df[keycols].astype(str).itertuples(index=False,name=None)))
    attempted=keys(p)
    numeric=keys(p.loc[p.real_numeric_ROI9])
    already_strict=keys(c.loc[c.strict_photo_gate])
    q["already_attempted_photo_cell"]=[k in attempted for k in
               q[keycols].astype(str).itertuples(index=False,name=None)]
    q["already_numeric_ROI9_cell"]=[k in numeric for k in
               q[keycols].astype(str).itertuples(index=False,name=None)]
    q["already_strict_photo_cell"]=[k in already_strict for k in
               q[keycols].astype(str).itertuples(index=False,name=None)]
    q["new_numeric_support_opportunity"]=~q.already_numeric_ROI9_cell
    eligible=q.loc[q.new_numeric_support_opportunity].copy()
    eligible=eligible.sort_values(["area_km2","land_id","year","cell_id"],
                                  ascending=[True,True,True,True])
    previous24k_upper=1321
    require(len(eligible)<=previous24k_upper,
            "whole photo pool cannot have more new cells than 24,576 subset")
    oldyears=keys(p)
    pair_keys=set(map(tuple,keep[["land_id","year"]].itertuples(index=False,name=None)))
    photo_years=set(map(tuple,p[["land_id","year"]].drop_duplicates().itertuples(index=False,name=None)))
    numeric_years=set(map(tuple,p.loc[p.real_numeric_ROI9,["land_id","year"]]
                         .drop_duplicates().itertuples(index=False,name=None)))
    rows={}
    for name,mask in [
        ("small_20km2_or_less",q.area_km2.le(20)),
        ("intermediate_20_to_100km2",q.area_km2.gt(20)&q.area_km2.le(100)),
        ("large_above100km2",q.area_km2.gt(100))
    ]:
        subset=q.loc[mask]
        rows[name]={"supported_islands":int(subset.land_id.nunique()),
                    "tier09_qualified_cells":int(len(subset)),
                    "full_pool_unmeasured_cells":int(subset.new_numeric_support_opportunity.sum()),
                    "full_pool_non_strict_cells":int((~subset.already_strict_photo_cell).sum())}
    result={
        "verified_tier09_light_islands":len(keep_ids),
        "verified_tier09_light_island_years":len(pair_keys),
        "previously_attempted_total_exact_photo_IDs":len(p),
        "previously_genuine_numeric_original_RGB_ROI9":int(p.real_numeric_ROI9.sum()),
        "previously_attempted_source_cells":len(attempted),
        "previously_numeric_source_cells":len(numeric),
        "previously_strict_source_cells":len(already_strict),
        "previously_strict_source_islands":int(c.loc[c.strict_photo_gate,"land_id"].nunique()),
        "previously_spatially_regressable_source_islands":133,
        "previously_spatially_regressable_source_siteyears":3315,
        "actual_overlapping_tier09_photo_IDs":len(duplicates),
        "tier09_qualified_existingVIIRS_cells":len(q),
        "tier09_qualified_cells_already_attempted":int(q.already_attempted_photo_cell.sum()),
        "tier09_qualified_cells_with_existing_true_numeric_ROI9":int(q.already_numeric_ROI9_cell.sum()),
        "tier09_qualified_cells_with_existing_strict_numeric_geography":int(q.already_strict_photo_cell.sum()),
        "tier09_qualified_cells_without_existing_true_numeric_ROI9":int(q.new_numeric_support_opportunity.sum()),
        "tier09_qualified_cells_not_already_strict":int((~q.already_strict_photo_cell).sum()),
        "tier09_islands_with_no_prior_numeric_ROI9":len(keep_ids-set(p.loc[p.real_numeric_ROI9,"land_id"])),
        "tier09_light_islandyears_with_no_prior_numeric_ROI9":len(pair_keys-numeric_years),
        "tier09_light_islandyears_with_no_prior_attempts":len(pair_keys-photo_years),
        "tier09_distinct_islands_with_unmeasured_qualified_cells":int(eligible.land_id.nunique()),
        "tier09_numeric_support_opportunity_by_true_area":rows,
    }
    return result,eligible


def run(photos,cells,fullreport,tier09,tier09report,pairs,out):
    out=Path(out)
    require(not out.exists(),"refuse to overwrite previous result")
    meta=json.loads(Path(fullreport).read_text())
    g=meta["deduplicated_geography_files"]
    require(g["photo_id_count"]==112428 and g["genuine_numeric_ROI9"]==40914 and
            sha(photos)==g["photo_id_SHA256"] and
            sha(cells)==g["cell_numeric_SHA256"],
            "current full RGB photo and original source-cell SHA mismatch")
    tier=json.loads(Path(tier09report).read_text())
    require(tier["schema"]=="ALAN_tier09_moth_sameyear_VIIRS_existingcell_source_gate_v1" and
            tier["source_roster_photo_count"]==78214 and
            sha(tier09)==tier["candidate_status_SHA256"] and
            tier["new_RGB_downloads"]==0 and tier["new_numeric_ROI9"]==0,
            "previous verified tier09 candidate/VIIRS receipt drift")
    dtype={"photo_id":str,"land_id":str,"year":str,"target_id":str,
           "cell_id":str,"species":str}
    p=pd.read_csv(photos,compression="gzip",dtype=dtype)
    c=pd.read_csv(cells,compression="gzip",dtype=dtype)
    n=pd.read_csv(tier09,compression="gzip",dtype=dtype)
    y=pd.read_csv(pairs,dtype={"land_id":str,"year":str})
    result,eligible=evaluate(p,c,n,y,meta)
    out.mkdir(parents=True)
    eligible.to_csv(out/"tier09_existingVIIRS_precolour_cells_not_yet_numeric_full112428.csv.gz",
                    index=False,compression="gzip")
    report={"schema":SCHEMA,"full_source_photo_ID_SHA256":sha(photos),
            "full_source_cell_SHA256":sha(cells),
            "tier09_original_photo_status_SHA256":sha(tier09),
            "additional_RGB_photo_attempts":0,"additional_numeric_ROI9":0,
            "not_original_photo_RGB_or_pollen_or_causal_effect":True,
            "no_old_numeric_observationID_new_UUID_crosswalk_claim":True,
            "no_GloBI_positive_means_exclusive_moth_or_nighttime_dependence":True,
            "supported_cells_not_new_independent_islands":True,
            "result":result}
    (out/"tier09_vs_entire_full112428_source_geography.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_TIER09_VS_ENTIRE_SIGNED_112428_GEOGRAPHY_VERIFIED",
          json.dumps(result,sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    for k in ("photos","cells","fullreport","tier09","tier09report","pairs","out"):
        a.add_argument("--"+k,required=True,type=Path)
    x=a.parse_args()
    run(x.photos,x.cells,x.fullreport,x.tier09,x.tier09report,x.pairs,x.out)
