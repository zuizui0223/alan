#!/usr/bin/env python3
"""Original independent 2015–2022 NPMS habitat-change declarations: schema & ID gate.

Source DOI 10.5285/a12b5c70-c0e2-48cf-b0b7-3e95791feca9.
The source is a distinct *published extract from the SAME NPMS monitoring
programme*, not an independently randomized or wholly separate field study.

Read only original byte-pinned ZIPs and report file schema, aggregate support
and exact record identifiers.  Never infer that 'no change recorded' means
ecologically unchanged; the response is *surveyor-reported habitat change*.
Never forward-fill missing habitat; never merge on rounded geographic positions.
No volunteer, plot ID, precise coordinates, free text or raw sample row exported.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

DOI="10.5285/a12b5c70-c0e2-48cf-b0b7-3e95791feca9"
SCHEMA="alan_npms_habitat_change_2015_2022_source_schema_gate_v1"
SOURCE_2025_SHA="eidxpinned_7files_2026-11-10"
SOURCE_MEMBERS={
    "data/occurrences_2015to2025.csv":(18724281,"616d29d2"),
    "data/npmshabitatlookup_fine2broad.csv":(1411,"60a53016"),
    "data/sampleattributes_2015to2025.csv":(12159752,"9c9465a3"),
    "data/npmshabitatlookup_surv2broad.csv":(1848,"4b2a1f88"),
    "data/dominscores.csv":(433,"64c7fb11"),
    "data/locations_2015to2025.csv":(200780,"d5d18b79"),
    "data/sampleinfowithlatlong_2015to2025.csv":(4095902,"5e52f384")
}
LIMIT=2000000

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as inp:
        for b in iter(lambda:inp.read(2<<20),b""):h.update(b)
    return h.hexdigest()

def safe_columns(df):
    # Header metadata is public schema, no rows or IDs are emitted.
    return [str(c) for c in df.columns]

def safe_values(col,limit=20):
    v=col.astype(str).str.strip()
    v=v[v.notna() & v.ne("") & v.ne("nan")]
    counts=v.value_counts()
    # Explicitly avoid outputting any potentially identifying/free-text
    # per-site values, IDs or coordinates.
    if (len(counts)>limit or counts.index.to_series().str.len().max()>45):
        return {"levels":int(len(counts)), "omitted":"not a small coded categorical field"}
    return {"levels":int(len(counts)),
            "value_counts":{str(k):int(n) for k,n in counts.items()}}

def is_coded_change(col):
    words=("changed","change","transition","different")
    negative=("exchange","interchange")
    name=str(col).lower()
    return any(x in name for x in words) and not any(x in name for x in negative)

def inspect(habitat_zip,npms_zip):
    result={
       "schema":SCHEMA,
       "habitat_dataset_doi":DOI,
       "independent_species_response_from_other_study":False,
       "same_NPMS_program_not_independent_external_replication":True,
       "habitat_source_zip_sha256":sha(habitat_zip),
       "newer_2015_2025_source_zip_sha256_not_used_as_fixed_id":sha(npms_zip),
       "source_original_habitat_rows_not_published":True,
       "source_plot_coordinates_not_published":True,
       "sample_habitat_effect_or_ALAN_moderator_fitted":False,
       "status":"source_only_preflight",
    }
    with zipfile.ZipFile(npms_zip) as old:
        assert old.testzip() is None
        for p,(size,crc) in SOURCE_MEMBERS.items():
            x=old.getinfo(p)
            assert x.file_size==size and f"{x.CRC:08x}"==crc,p
        with old.open("data/sampleinfowithlatlong_2015to2025.csv") as f:
            source=pd.read_csv(f,dtype=str,keep_default_na=False)
    assert len(source)==27611 and source.id.is_unique
    result["old_NPMS_source_sample_rows"]=len(source)
    result["old_NPMS_true_location_count"]=source.location_id.nunique()
    with zipfile.ZipFile(habitat_zip) as z:
        assert z.testzip() is None
        members=[x for x in z.infolist() if not x.is_dir()]
        result["source_files"]=[]
        for member in members:
            entry={"name":member.filename,
                   "size":member.file_size,
                   "original_zip_CRC32":f"{member.CRC:08x}"}
            if member.filename.lower().endswith((".csv",".tsv")):
                raw=z.read(member)
                sep="\t" if member.filename.lower().endswith(".tsv") else ","
                table=pd.read_csv(io.BytesIO(raw),dtype=str,
                                  sep=sep,keep_default_na=False,
                                  low_memory=False,encoding_errors="replace")
                entry.update(columns=safe_columns(table),
                             source_row_count=len(table),
                             source_unique_row_count=len(table.drop_duplicates()))
                exact_identifiers={}
                for candidate,other in (
                    ("location_id","location_id"),("plot_id","location_id"),
                    ("sample_id","id"),("survey_id","id"),("id","id")):
                    if candidate in table.columns:
                        nonempty=table[candidate].loc[table[candidate].ne("")]
                        truth=set(source[other])-set([""])
                        exact_identifiers[candidate]={
                            "original_NPMS_compared_field":other,
                            "rows_with_exact_NPMS_id":int(nonempty.isin(truth).sum()),
                            "distinct_source_keys_with_NPMS_match":int(nonempty[nonempty.isin(truth)].nunique()),
                            "distinct_nonempty_habitat_keys":int(nonempty.nunique())
                        }
                entry["exact_NPMS_record_identifier_overlap"]=exact_identifiers
                entry["candidate_change_flags"]={}
                for c in table.columns:
                    if is_coded_change(c):
                        entry["candidate_change_flags"][c]=safe_values(table[c])
                datecols=[c for c in table.columns if "date" in c.lower() or
                          "year" in c.lower()]
                entry["date_columns"]=datecols
                result["source_files"].append(entry)
            else:
                result["source_files"].append(entry)
    result["potentially_joinable_files"]=[
        r["name"] for r in result["source_files"] if
        any(q["rows_with_exact_NPMS_id"]>0 for q in
            r.get("exact_NPMS_record_identifier_overlap",{}).values())
    ]
    result["status"]=(
        "source_original_habitat_schema_and_IDs_verified"
        if result["potentially_joinable_files"] else
        "no_exact_key_overlap_STOP"
    )
    result["next_model_stopped_until"]=[
        "Source verifies whether any unique, explicit yes/no change field is observed",
        "Source IDs can be joined by true plot or original sample key with full accuracy",
        "When no change flag was filed, missingness stays UNKNOWN, never zero",
        "Original 2015-2022 habitat records must be matched to 2025-source sample IDs after version drift audit",
        "Species-level positive cover and numerical PCNL present only in old fixed source and must never be overwritten"
    ]
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--habitat-2022",required=True)
    p.add_argument("--npms-2025",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    result=inspect(args.habitat_2022,args.npms_2025)
    dest=Path(args.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("ALAN_NPMS_OFFICIAL_HABITAT_CHANGE_SOURCE",json.dumps(result,
           ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
