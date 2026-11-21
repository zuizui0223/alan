#!/usr/bin/env python3
"""NPMS 2015–25 actual source-survey habitats and explicit 'no species' flag gate.

The documented species-level MOIST signal is not evidence that the surveyed
sites are literally wetlands. This audit tests that distinction with original
volunteer NPMS habitat labels and a fixed repeated-cover source cohort. It
never labels species absent because a list lacks that species.
No volunteer identifiers or exact location IDs are published.
"""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

MEMBERS={
    "data/occurrences_2015to2025.csv":(18724281,"616d29d2"),
    "data/npmshabitatlookup_fine2broad.csv":(1411,"60a53016"),
    "data/sampleattributes_2015to2025.csv":(12159752,"9c9465a3"),
    "data/npmshabitatlookup_surv2broad.csv":(1848,"4b2a1f88"),
    "data/dominscores.csv":(433,"64c7fb11"),
    "data/locations_2015to2025.csv":(200780,"d5d18b79"),
    "data/sampleinfowithlatlong_2015to2025.csv":(4095902,"5e52f384")
}
SCHEMA="alan_npms_original_habitat_and_recording_status_support_v1"
def guard(x,msg):
    if not bool(x):raise RuntimeError("NPMS_HABITAT_SOURCE_STOP: "+msg)
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as s:
        for chunk in iter(lambda:s.read(2<<20),b""):h.update(chunk)
    return h.hexdigest()
def top(s):
    return {str(k):int(v) for k,v in s.value_counts(dropna=False).items()}
def read(z):
    with zipfile.ZipFile(z) as zp:
        guard(zp.testzip() is None,"published ZIP CRC failed")
        for p,(size,crc) in MEMBERS.items():
            i=zp.getinfo(p)
            guard(i.file_size==size and f"{i.CRC:08x}"==crc,"source member CRC/size mismatch: "+p)
        out={}
        for n,p in {
            "sample":"data/sampleinfowithlatlong_2015to2025.csv",
            "attrs":"data/sampleattributes_2015to2025.csv",
            "lookup":"data/npmshabitatlookup_surv2broad.csv",
            "fine_lookup":"data/npmshabitatlookup_fine2broad.csv",
            "occ":"data/occurrences_2015to2025.csv"
        }.items():
            with zp.open(p) as f:out[n]=pd.read_csv(f,dtype=str,keep_default_na=False)
    guard(len(out["sample"])==27611 and len(out["attrs"])==293561 and
          len(out["occ"])==271701,"NPMS row counts differ from locked original")
    return out

def analyze(z,old,receipt,return_matched=False):
    src=read(z)
    sample=src["sample"]
    att=src["attrs"]
    occ=src["occ"]
    look=src["lookup"]
    assert {"id","survey_id","location_id","date_start"}.issubset(sample.columns)
    assert {"sample_id","caption","text_value"}.issubset(att.columns)
    assert {"NPMS_surveyed","NPMS_broad_habitat"}.issubset(look.columns)
    ids=set(sample.id)
    guard(sample.id.is_unique,"original duplicated sample key")
    a=att.loc[att.sample_id.isin(ids)].copy()
    hab=a[a.caption=="NPMS Habitat"].copy()
    nospecies=a[a.caption=="No species observed"].copy()
    hab["raw"]=hab.text_value.astype(str).str.strip()
    nospecies["raw"]=nospecies.text_value.astype(str).str.strip()
    # Preserve *original* survey habitat labels and supplied broad mappings.
    look["raw"]=look.NPMS_surveyed.astype(str).str.strip()
    mappings=look.groupby("raw").NPMS_broad_habitat.agg(lambda x:
                sorted(set(z for z in x if z))).to_dict()
    def map_exact(s):
        opts=mappings.get(s,[])
        return opts[0] if len(opts)==1 else ""
    hab["broad"]=hab.raw.map(map_exact)
    # If multiple habitats are recorded within one sample, never arbitrarily
    # select one value to rescue sample size.
    by_sample=hab.groupby("sample_id",sort=False).agg(
        n_labels=("raw","nunique"),
        n_broad=("broad",lambda x:len(set(v for v in x if v))),
        broad=("broad",lambda x:next(iter(set(v for v in x if v)),"")),
        raw=("raw",lambda x:next(iter(set(v for v in x if v)),"")))
    by_sample.loc[(by_sample.n_labels!=1)|(by_sample.n_broad>1),"broad"]=""
    sm=sample[["id","survey_id","location_id","date_start"]].copy()
    sm["date"]=pd.to_datetime(sm.date_start,utc=True,errors="coerce")
    sm["year"]=sm.date.dt.year
    sm["season"]=np.select([sm.date.dt.month.isin([4,5,6]),
                            sm.date.dt.month.isin([7,8,9])],
                           ["spring","summer"],default="other")
    sm["method"]=sm.survey_id.str.strip()
    sm=sm.join(by_sample[["broad","raw"]],on="id")
    sm=sm.loc[sm.location_id.ne("")&sm.year.between(2015,2025)]
    sy=sm[["location_id","year","season","method","broad","raw"]]
    group=sy.groupby(["location_id","year","season","method"],sort=False,dropna=False).agg(
        n_survey_records=("broad","size"),
        n_nonempty_habitat=("broad",lambda x:sum(str(v).strip()!="" and pd.notna(v) for v in x)),
        original_classes=("broad",lambda x:sorted(set(str(v) for v in x if pd.notna(v) and str(v).strip()))),
        original_raw_classes=("raw",lambda x:sorted(set(str(v) for v in x if pd.notna(v) and str(v).strip()))),
    ).reset_index()
    group["valid_unique_habitat"]=group.original_classes.map(len).eq(1)
    group["broad"]=group.original_classes.map(lambda x:x[0] if len(x)==1 else "")
    # Join only matching observed year, method and broad season; some NPMS
    # plot-year strata have >1 visit, and conflicting habitats stay unresolved.
    original_receipt=json.loads(Path(receipt).read_text())
    guard(original_receipt.get("source_original_positive_cover_table_SHA256")==sha(old),
          "source positive-cover gzip SHA differs from exact author-validated audit")
    orig=pd.read_csv(old,dtype={"location_id":str,"survey_method_id":str,"season":str,"species":str},low_memory=False)
    # The source-level CSV is broader than the 61,555 source-positive original
    # outcome rows. Reapply the SAME frozen source panel eligibility used in
    # the old original 2015-25 positive-cover fixed-effect study.
    req={"location_id","survey_method_id","season","species","year",
         "original_annual_priorPCNL","source_domin_encodings"}
    guard(req.issubset(orig.columns),"missing original published panel columns")
    keys=["location_id","survey_method_id","season","species"]
    orig["panel_key"]=orig[keys].agg("|".join,axis=1)
    g=orig.groupby("panel_key",sort=False).agg(
        yrs=("year","nunique"),levels=("original_annual_priorPCNL","nunique"),
        min_light=("original_annual_priorPCNL","min"),
        max_light=("original_annual_priorPCNL","max"),
        domin=("source_domin_encodings","nunique"))
    keep=g.loc[(g.yrs>=2)&(g.levels>=2)&
               (g.max_light-g.min_light>=.02)&(g.domin==1)].index
    orig=orig.loc[orig.panel_key.isin(keep)].copy()
    guard(len(orig)==61555,"old original fixed sample cohort changed after locked panel gate")
    guard(orig.location_id.nunique()==1138 and orig.species.nunique()==686,
          "old NPMS original plot/taxon support changed")
    guard(orig.year.between(2015,2025).all(),"invalid original year")
    joined=orig.merge(group[["location_id","year","season","method","valid_unique_habitat","broad"]],
        left_on=["location_id","year","season","survey_method_id"],
        right_on=["location_id","year","season","method"],
        how="left",validate="many_to_one")
    guard(len(joined)==61555,"source duplicated after habitat join")
    joined["valid_unique_habitat"]=joined.valid_unique_habitat.fillna(False)
    good=joined.loc[joined.valid_unique_habitat].copy()
    # By original source-habitat category, quantify geographical and repeated
    # positive-cover overlap. This is NOT a habitat-specific effect estimate.
    h=good.groupby("broad",dropna=False).agg(rows=("species","size"),
        plots=("location_id","nunique"),species=("species","nunique"),
        years=("year","nunique")).reset_index()
    study={"schema":SCHEMA,"published_source_doi":"10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8",
        "source_CRC_verified":True,"source_zip_outer_sha_unstable_but_member_CRCS_fixed":True,
        "source_cover_sha256":sha(old),"original_cover_rows":61555,
        "attr_caption_counts":top(att.caption),
        "npms_original_habitat_raw_top":top(hab.raw) if len(hab) else {},
        "npms_original_habitat_mapped_to_broad_n":int(hab.broad.ne("").sum()),
        "npms_original_habitat_n_samples":int(hab.sample_id.nunique()),
        "npms_original_habitat_broad_top":top(hab.broad),
        "broad_lookup_labels":list(sorted(set(look.NPMS_broad_habitat))),
        "sample_no_species_flag_unique_values":top(nospecies.raw),
        "sample_no_species_flag_source_records":len(nospecies),
        "source_occurrence_status_codes":top(occ.record_status),
        "source_cover_rows_with_unique_recorded_broad_habitat":len(good),
        "source_cover_unique_broad_habitat_plot_count":int(good.location_id.nunique()),
        "positive_cover_broad_classes":
            [{k:(int(row[k]) if k in ("rows","plots","species","years") else str(row[k]))
              for k in h.columns} for _,row in h.iterrows()],
        "target":"independent habitat × MOIST × historical prior PCNL recorded-cover support",
        "plant_local_absence_fitted":False,
        "habitat_specific_effect_fitted":False,
        "reason_no_zeros":"A missing species row is not a species-specific absence. 'No species observed' is a survey-level attribute until protocol-specific negative ascertainment is verified.",
        "status":"source_habitat_coverage_audited_no_effect"}
    return (study,good) if return_matched else study

def cli():
    p=argparse.ArgumentParser()
    p.add_argument("--npms",required=True)
    p.add_argument("--cover",required=True)
    p.add_argument("--receipt",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    out=analyze(a.npms,a.cover,a.receipt)
    dest=Path(a.out);dest.parent.mkdir(exist_ok=True,parents=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("NPMS_OFFICIAL_HABITAT_SOURCE_GATE",json.dumps({
        "study_status":out["status"],
        "observed_broad_habitats":out["positive_cover_broad_classes"],
        "no_species_flag_values":out["sample_no_species_flag_unique_values"],
        "mapped_habitat_source_samples":out["npms_original_habitat_mapped_to_broad_n"],
        "source_record_status":out["source_occurrence_status_codes"]
    },ensure_ascii=False,sort_keys=True))

if __name__=="__main__":cli()
