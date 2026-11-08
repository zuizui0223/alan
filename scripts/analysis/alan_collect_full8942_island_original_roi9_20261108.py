#!/usr/bin/env python3
"""SHA-authenticated full 8,942 true NEW island photo ROI9 acquisition audit.

Read only previous 5×256 signed original photos and 15 new warm-group original
ROI outcomes. Require EXACT photo ID and original observation UUID equality
with the fixed 8,942 source roster. Technical ROI failures remain unknown.
Contrast what source-only 5-photo/3-WCVP-species threshold *could* support
with what actually meets 5 attempts, 3 numeric ROI9 successes from 3 distinct
photographed accepted taxa. Compare island-breadth and separately frozen
site-year-focused 512 acquisition subset on real pixels WITHOUT reinference.
This is not the original ALAN historical-main ecological regression.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np,pandas as pd
from scripts.acquisition import alan_bulk_7662_original_roi9_source_20261108 as source

FIRST_EXPECTED={1:92,2:91,3:98,4:106,5:105}
FIRST_SOURCE_SHA="7789c5aa1916e5b1ff915d217744a80b990795cb24e362cff20007dda96ee94b"
DENSE512_SHA="172680e798ba9711ee149fd75fc84601e55137c17a818c028aa36ea24c309697"
MODEL_BLOB=source.MODEL_BLOB
N=8942

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN8942_ROI9_SOURCE_AUDIT: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def load_original_file(file,expected,n_expected,sig_required=True):
    p=Path(file)
    df=pd.read_csv(p,compression="gzip",dtype=str,keep_default_na=False)
    require(len(df)==n_expected and df.photo_id.is_unique and df.observation_id.is_unique,
            "Previously signed original ROI output lost independent observation UUID")
    require("strict9" in df and "strict9_white" in df,
            "Original FCP numerical nine-palette source fields missing")
    df["real_numeric_ROI9"]=df.strict9.str.lower().eq("true")
    require(int(df.real_numeric_ROI9.sum())==expected,
            "Previously signed original nine-palette colour success count drift")
    vals=pd.to_numeric(df.loc[df.real_numeric_ROI9,"strict9_white"],errors="coerce")
    require(vals.between(0,1).all(),
            "Classified original photographic ROI9 white fraction invalid")
    require(df.loc[~df.real_numeric_ROI9,"strict9_white"].eq("").all(),
            "Unclassified colour was silently filled")
    return df[["photo_id","observation_id","real_numeric_ROI9"]]

def load_prior(root,n):
    root=Path(root)
    if n==1:
        file=root/"finished"/"new256_independent_island_actual_original_ROI9.csv.gz"
        rec=root/"finished"/"ALAN_new256_original_photo_ROI9_source_receipt.json"
        require(sha(file)==FIRST_SOURCE_SHA,"Original signed first256 image pixel outcome changed")
        r=json.loads(rec.read_text())
        require(r.get("schema")=="ALAN_new256_original_island_unique_UUID_real_FCP_ROI9_v1"
                and r["original_photo_colour_output_SHA256"]==FIRST_SOURCE_SHA and
                r["original_nine_fraction_strict_colour_success"]==FIRST_EXPECTED[1],
                "First signed original source receipt mismatch")
    else:
        file=root/"finished"/f"ALAN_original_novel_island_photo_ROI9_batch_{n:02d}.csv.gz"
        rec=root/"finished"/"ALAN_original_novel_island_ROI9_batch_receipt.json"
        r=json.loads(rec.read_text())
        require(r.get("schema")=="ALAN_novel_source_island_original256_FCP_ROI9_batch_v1"
                and r["source_batch"]==n and
                r["numeric_original_nine_palette_success"]==FIRST_EXPECTED[n] and
                r["measured_original_colours_sha256"]==sha(file),
                "Original prior batch photo source SHA/ROI9 count mismatch")
    df=load_original_file(file,FIRST_EXPECTED[n],256)
    df["original_source_batch"]=n
    return df,{"source_batch":n,"genuine_original_photo_attempts":len(df),
               "numerical_original_ROI9_success":FIRST_EXPECTED[n],
               "original_photo_output_SHA256":sha(file)}

def load_bulk(root,g):
    root=Path(root)
    file=root/"finished"/f"real_new_original_island_ROI9_group{g:02d}.csv.gz"
    rec=root/"finished"/"real_original_island_bulk512_group_receipt.json"
    r=json.loads(rec.read_text())
    require(r["group"]==g and
            r["original_source_photo_attempts"]==source.group_size(g) and
            r["new_original_independent_observation_UUIDs"]==source.group_size(g) and
            r["original_numeric_colour_file_SHA256"]==sha(file) and
            r["source_method_original_git_blob"]==MODEL_BLOB and
            r["source_image_bytes_persisted"] is False and
            r["unknown_photo_colour_was_not_imputed"] is True,
            "Missing or altered 15-group actual FCP ROI9 original-image receipt")
    frame=load_original_file(file,r["original_nine_palette_ROI9_true_success"],
                             source.group_size(g))
    return frame,{"warm_group":g,
                  "source_batch_pair":list(source.group_batches(g)),
                  "genuine_original_photo_attempts":len(frame),
                  "numerical_original_ROI9_success":int(frame.real_numeric_ROI9.sum()),
                  "original_photo_output_SHA256":sha(file)}

def count_source_sites(df):
    require(df.observation_id.is_unique and df.photo_id.is_unique,
            "Actual original photo pixel dataset pseudo-replicates observation UUID")
    attempted=df.groupby("target_id",observed=True).agg(
       n_photo_observation_UUIDs=("observation_id","nunique"),
       observed_source_accepted_taxa=("wcvp_accepted_plant_name_id","nunique"),
       actually_measured_ROI9=("real_numeric_ROI9","sum"))
    measured=df.loc[df.real_numeric_ROI9]
    accepted=measured.groupby("target_id",observed=True).wcvp_accepted_plant_name_id.nunique()
    attempted["measured_distinct_accepted_taxa"]=accepted.reindex(
         attempted.index,fill_value=0).astype(int)
    eligible_source=(attempted.n_photo_observation_UUIDs.ge(5)&
                     attempted.observed_source_accepted_taxa.ge(3))
    eligible_actual=(eligible_source&attempted.actually_measured_ROI9.ge(3)&
                     attempted.measured_distinct_accepted_taxa.ge(3))
    return {
      "original_photo_observation_UUIDs":len(df),
      "actual_nine_palette_measured_photos":int(df.real_numeric_ROI9.sum()),
      "actual_original_ROI9_unknown_photos":int((~df.real_numeric_ROI9).sum()),
      "distinct_GSHHG_islands":int(df.land_id.nunique()),
      "source_photo_cell_years":int(df.target_id.nunique()),
      "source_opportunity_site_years_5attempts_3taxa":int(eligible_source.sum()),
      "actual_5attempt_3numericROI9_3observed_accepted_taxa_site_years":int(eligible_actual.sum()),
      "site_years_with_1plus_numericROI9":int(
           attempted.actually_measured_ROI9.ge(1).sum()),
    }

def run(schedule,source_receipt,first,older,bulk,dense512,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting SHA-frozen completed original photo audit")
    master,_=source.source_read(schedule,source_receipt)
    require(len(older)==4 and set(older)=={2,3,4,5} and
            len(bulk)==15 and set(bulk)==set(range(15)),
            "Incomplete original authenticated first 5 or last 15 warm worker sources")
    frames=[];proof=[]
    for n,root in [(1,first),*sorted(older.items())]:
        d,r=load_prior(root,n)
        frames.append(d);proof.append(r)
    for group,root in sorted(bulk.items()):
        d,r=load_bulk(root,group)
        frames.append(d);proof.append(r)
    measured=pd.concat(frames,ignore_index=True)
    require(len(measured)==N and measured.photo_id.is_unique and
            measured.observation_id.is_unique,
            "Original 8942 new source photos lost/duplicated across 15+5 workers")
    expected=master[["photo_id","observation_id","source_batch","land_id",
                     "target_id","wcvp_accepted_plant_name_id"]]
    joined=expected.merge(measured,on=["photo_id","observation_id"],
                          validate="one_to_one",how="inner")
    require(len(joined)==N and set(joined.photo_id)==set(master.photo_id) and
            set(joined.observation_id)==set(master.observation_id),
            "Original FCP ROI results do not correspond to exact source photo/UUID roster")
    all_stats=count_source_sites(joined)
    require(all_stats["original_photo_observation_UUIDs"]==N and
            all_stats["source_opportunity_site_years_5attempts_3taxa"]==303 and
            all_stats["source_photo_cell_years"]==4776 and
            all_stats["distinct_GSHHG_islands"]==369,
            "Original 8942 source photo/taxon/sites changed")
    first5=joined.loc[joined.source_batch.astype(int).le(5)]
    first7=joined.loc[joined.source_batch.astype(int).le(7)]
    first5_stats=count_source_sites(first5)
    first7_stats=count_source_sites(first7)
    require(first5_stats["actual_nine_palette_measured_photos"]==492 and
            first5_stats["actual_5attempt_3numericROI9_3observed_accepted_taxa_site_years"]==1 and
            first5_stats["source_opportunity_site_years_5attempts_3taxa"]==3 and
            first7_stats["source_opportunity_site_years_5attempts_3taxa"]==6,
            "First 1280 actual ROI and first1792 source photo site counts changed")
    require(sha(dense512)==DENSE512_SHA,
            "Original separately frozen site-density512 source photo IDs changed")
    density=pd.read_csv(dense512,compression="gzip",dtype=str,keep_default_na=False)
    require(len(density)==512 and density.photo_id.is_unique and
            density.observation_id.is_unique,
            "Original site-density source photos/UUIDs duplicated")
    dense_set=set(density.observation_id)
    require(dense_set.isdisjoint(set(first5.observation_id)) and
            dense_set.issubset(set(joined.observation_id)),
            "Original site-density wave observation UUIDs leak earlier pixels")
    dense=joined.loc[joined.observation_id.isin(set(first5.observation_id)|dense_set)]
    require(len(dense)==1792,"Actual source colour comparison not same 1792 observations")
    dense_stats=count_source_sites(dense)
    require(dense_stats["source_opportunity_site_years_5attempts_3taxa"]==60,
            "Pre-colour 60 strict source-site opportunity gain changed")
    support={
       "first1280_island_breadth":first5_stats,
       "first1792_breadth_source":first7_stats,
       "first1280_plus512_island_site_densification":dense_stats,
       "all8942_completed_original_island_source":all_stats
    }
    counts=joined.groupby("source_batch",observed=True).agg(
        original_photo_observation_UUIDs=("observation_id","nunique"),
        numeric_original_ROI9=("real_numeric_ROI9","sum"),
        independent_GSHHG_islands=("land_id","nunique"),
        source_site_years=("target_id","nunique")).reset_index()
    out.mkdir(parents=True)
    outfile=out/"ALAN_real_original8942_ROI9_colour_source_batch_counts.csv"
    counts.to_csv(outfile,index=False)
    receipt={
      "schema":"ALAN_real_full8942_original_island_ROI9_source_SHA_merge_v1",
      "source_original_full_8942_photo_schedule_SHA256":source.SOURCE_SHA,
      "new_original_photo_IDs_and_independent_UUIDs_verified":N,
      "new_original_nine_palette_true_success":all_stats["actual_nine_palette_measured_photos"],
      "original_colour_unclassifiable_photos":all_stats["actual_original_ROI9_unknown_photos"],
      "original_source_islands":369,
      "original_source_photo_cell_years":4776,
      "site_year_sample_support":support,
      "per_source_batch_summary_SHA256":sha(outfile),
      "original_FCP_evaluator_blob_SHA256":MODEL_BLOB,
      "number_original_SHA_signed_image_results_included":len(proof),
      "per_original_image_source_provenance":proof,
      "original_photo_image_bytes_reused_without_new_download_or_inference":True,
      "unknown_original_ROI9_never_inferred_as_white":True,
      "new_source_island_targeted_NOT_unbiased_global_flora_sample":True,
      "original_2019_2021_2023_2024_2025_historical_main_ALAN_model_unchanged":True,
    }
    (out/"ALAN_real_8942_original_island_photo_ROI9_final_receipt.json").write_text(
       json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_8942_ORIGINAL_ISLAND_NUMERICAL_ROI9_AUDITED",
          json.dumps({k:receipt[k] for k in (
           "new_original_photo_IDs_and_independent_UUIDs_verified",
           "new_original_nine_palette_true_success",
           "original_colour_unclassifiable_photos","site_year_sample_support")
           },sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schedule",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--first",type=Path,required=True)
    for n in range(2,6):p.add_argument(f"--prior{n}",type=Path,required=True)
    for n in range(15):p.add_argument(f"--group{n}",type=Path,required=True)
    p.add_argument("--densification512",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    q=p.parse_args()
    run(q.schedule,q.source_receipt,q.first,
        {n:getattr(q,f"prior{n}") for n in range(2,6)},
        {n:getattr(q,f"group{n}") for n in range(15)},
        q.densification512,q.out)
