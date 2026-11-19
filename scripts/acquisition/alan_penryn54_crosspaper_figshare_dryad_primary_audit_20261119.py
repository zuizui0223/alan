#!/usr/bin/env python3
"""Original dataset source audit: one 54-mesocosm programme, three papers.

Programs:
Bennie et al. 2015 resource-mediated aphids; Bennie et al. 2018 predator
presence × light; Anic et al. 2022 plant biomass and shoot:root allocation.
Physical lights: 3 spectra × 3 trophic compositions × 6 mesocosms = 54.
This audit tries official Dryad and Figshare APIs, downloads bounded public
original files when possible, publishes file SHA and anonymised *schema only*,
never original unit rows, treatment allocations or plot identifiers.

A green pipeline with SOURCE_STOP does NOT imply original rows were retrieved.
"""
import argparse
import csv
import hashlib
import io
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

FIGSHARE_ARTICLE=20304039
FIGSHARE_DOI="10.6084/m9.figshare.20304039"
DRYAD_DOI="10.5061/dryad.8217n76"
DRYAD_API="https://datadryad.org/api/v2"
FIGSHARE_API="https://api.figshare.com/v2"
# User is only requesting a public ecological study. Still prohibit exporting
# exact plot/site identifiers and individual biological observation rows.
MAX_BYTES=18_000_000
CANDIDATE=re.compile(r"(plot|site|meso|trough|unit|rep|block|year|season|light|troph|pred|herb|treat|spec)",re.I)
UNIT_CANDIDATE=re.compile(r"(plot|meso|trough|site|unit|block)",re.I)


def pull(url,limit=MAX_BYTES):
    request=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-open-biodiversity-source-reproducibility/2026-11-19",
        "Accept":"application/json,text/plain,text/csv,application/zip,application/octet-stream,*/*",
    })
    with urllib.request.urlopen(request,timeout=30) as response:
        data=response.read(limit+1)
        if len(data)>limit:
            raise ValueError("Original source exceeds bounded preflight")
        return data


def classify_original_table(header,rows):
    cols=[str(x).strip() if x is not None else "" for x in header]
    # Original header text and row counts are sufficient for source access
    # audit; actual contents of an identified mesocosm_ID column stay private.
    n=0
    distinct={str(i):set() for i,c in enumerate(cols) if CANDIDATE.search(c)}
    for row in rows:
        if n>=250_000:raise ValueError("Original source rows exceed preflight")
        if not any(x is not None and str(x).strip() for x in row):
            continue
        n+=1
        for s in distinct:
            i=int(s)
            if i<len(row):
                value=row[i]
                if value is not None and str(value).strip():
                    distinct[s].add(str(value).strip())
    return {
        "column_names_from_exact_original_file":cols,
        "nonempty_data_rows":n,
        "candidate_source_fields":[
            {"source_column_index":int(i),"source_column_name":cols[int(i)],
             "n_distinct_original_values_only":len(v)}
            for i,v in sorted(distinct.items(),key=lambda z:int(z[0]))
        ],
        "exact_original_plot_identifiers_exported":False,
        "original_observation_values_exported":False,
        "original_treatment_mean_or_SE_fitted":False,
    }


def csv_schema(original,delimiter=","):
    txt=original.decode("utf-8-sig",errors="replace")
    z=list(csv.reader(io.StringIO(txt),delimiter=delimiter))
    if len(z)<2:
        raise ValueError("Original source table has no data rows")
    return classify_original_table(z[0],z[1:])


def xlsx_schema(original):
    from openpyxl import load_workbook
    if not original.startswith(b"PK\x03\x04"):
        raise ValueError("Not original XLSX ZIP signature")
    wb=load_workbook(io.BytesIO(original),read_only=True,data_only=True)
    out=[]
    for ws in wb.worksheets:
        nonempty=[]
        for i, row in enumerate(ws.iter_rows(values_only=True),start=1):
            if not any(x is not None and str(x).strip() for x in row):
                continue
            nonempty.append((i,list(row)))
            if len(nonempty)>250_002:
                raise ValueError("Source workbook exceeds bounded row count")
        if len(nonempty)<2:
            out.append({"sheet_name":ws.title,"nonempty_rows":len(nonempty),
                        "original_source_table_unidentified":True})
            continue
        # Identify header among first 12 nonempty rows; prefer the first
        # plausible multiple-column source header, without guessing identities.
        idx=next((i for i,(_,r) in enumerate(nonempty[:12])
                  if sum(isinstance(c,str) and bool(c.strip()) for c in r)>=3),None)
        if idx is None:
            out.append({"sheet_name":ws.title,"nonempty_rows":len(nonempty),
                        "original_source_table_unidentified":True})
            continue
        d=classify_original_table(nonempty[idx][1],(r for _,r in nonempty[idx+1:]))
        d["sheet_name"]=ws.title
        d["source_header_original_excel_row"]=nonempty[idx][0]
        out.append(d)
    wb.close()
    return out


def safe_schema(name,b):
    suffix=Path(name).suffix.lower()
    if suffix==".csv":
        return {"kind":"csv","sheets":[csv_schema(b)]}
    if suffix==".tsv":
        return {"kind":"tsv","sheets":[csv_schema(b,delimiter="\t")]}
    if suffix==".xlsx":
        return {"kind":"xlsx","sheets":xlsx_schema(b)}
    if suffix==".zip" and b.startswith(b"PK\x03\x04"):
        # Bound total uncompressed size and number of members against ZIP bombs.
        sheets=[]
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            members=z.infolist()
            if len(members)>100:
                raise ValueError("Original ZIP has too many members")
            if sum(m.file_size for m in members)>30_000_000:
                raise ValueError("Original ZIP too large uncompressed")
            for f in members:
                if f.is_dir():continue
                if f.filename.lower().endswith((".csv",".tsv",".xlsx")):
                    sub=z.read(f)
                    d=safe_schema(f.filename,sub)
                    sheets.append({"archive_member":f.filename,
                                   "member_SHA256":hashlib.sha256(sub).hexdigest(),
                                   **d})
        return {"kind":"zip","members":sheets}
    return {"kind":"other","original_data_type_unparsed":suffix or "unknown"}


def get_figshare():
    out={
        "dataset_DOI":FIGSHARE_DOI,"original_article":FIGSHARE_ARTICLE,
        "metadata":"SOURCE_STOP","files":[],
        "source_authenticated_original_table_count":0,
        "original_light_meso_ids_jointly_matched_other_source":False,
    }
    try:
        metadata=pull(f"{FIGSHARE_API}/articles/{FIGSHARE_ARTICLE}",limit=2_000_000)
        doc=json.loads(metadata)
        observed=doc.get("doi","")
        if observed and FIGSHARE_DOI not in observed.lower():
            raise ValueError("Different Figshare DOI: refusing unexpected source")
        out["metadata"]="VERIFIED_PUBLISHED_FIGSHARE_API"
        out["title"]=doc.get("title")
        out["metadata_sha256"]=hashlib.sha256(metadata).hexdigest()
        for f in doc.get("files",[])[:50]:
            info={
                "original_file_id":f.get("id"),"original_file_name":f.get("name"),
                "declared_size":f.get("size"),
                "upstream_checksum_MD5":f.get("computed_md5") or f.get("supplied_md5"),
                "download_status":"SOURCE_STOP",
            }
            url=f.get("download_url","")
            size=f.get("size")
            if isinstance(url,str) and url.startswith("https://") and isinstance(size,int) and 0<size<=MAX_BYTES:
                try:
                    raw=pull(url)
                    info["download_status"]="VERIFIED_RAW_BYTES"
                    info["original_file_SHA256"]=hashlib.sha256(raw).hexdigest()
                    info["original_file_size_bytes"]=len(raw)
                    info["original_file_schema_only"]=safe_schema(f.get("name",""),raw)
                    if info["original_file_schema_only"]["kind"] in ("csv","tsv","xlsx","zip"):
                        out["source_authenticated_original_table_count"]+=1
                except (ValueError,ImportError,zipfile.BadZipFile,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
                    info["error_type"]=type(exc).__name__
                    if isinstance(exc,urllib.error.HTTPError):
                        info["http_status"]=exc.code
            out["files"].append(info)
    except (ValueError,KeyError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
        out["metadata_failure"]={"error_type":type(exc).__name__,
                                 **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})}
    return out


def get_dryad():
    out={"dataset_DOI":DRYAD_DOI,"metadata":"SOURCE_STOP",
         "files":[],"source_authenticated_original_table_count":0,
         "original_light_meso_ids_jointly_matched_other_source":False}
    doi=urllib.parse.quote("doi:"+DRYAD_DOI,safe="")
    for label,url in [
        ("dataset",f"{DRYAD_API}/datasets/{doi}"),
        ("versions",f"{DRYAD_API}/datasets/{doi}/versions"),
    ]:
        try:
            b=pull(url,limit=2_000_000)
            data=json.loads(b)
            out[label+"_metadata_status"]="VERIFIED_DRYAD_PUBLIC_API"
            out[label+"_metadata_sha256"]=hashlib.sha256(b).hexdigest()
            if label=="dataset":
                out["dataset_title"]=data.get("title","")
                out["dataset_version_metadata_keys"]=sorted(data.keys())
            if label=="versions":
                out["versions_metadata_keys"]=sorted(data.keys()) if isinstance(data,dict) else []
                embedded=data.get("_embedded",{}) if isinstance(data,dict) else {}
                ids=[]
                for v in embedded.values():
                    if isinstance(v,list):
                        for rec in v:
                            if isinstance(rec,dict) and isinstance(rec.get("id"),int):
                                ids.append(rec["id"])
                out["resolved_version_IDs"]=ids[:5]
        except (ValueError,KeyError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            out[label+"_metadata_status"]="SOURCE_STOP"
            out[label+"_error_type"]=type(exc).__name__
            if isinstance(exc,urllib.error.HTTPError):
                out[label+"_http_status"]=exc.code
    for ver in out.get("resolved_version_IDs",[])[:4]:
        try:
            blob=pull(f"{DRYAD_API}/versions/{ver}/files",limit=2_000_000)
            doc=json.loads(blob)
            out.setdefault("verified_versions_file_index",[]).append({
                "version_ID":ver,
                "metadata_SHA256":hashlib.sha256(blob).hexdigest(),
                "original_file_list_keys":sorted(doc.keys()) if isinstance(doc,dict) else [],
            })
            embedded=doc.get("_embedded",{}) if isinstance(doc,dict) else {}
            recs=[]
            for seq in embedded.values():
                if isinstance(seq,list):
                    recs.extend(x for x in seq if isinstance(x,dict))
            for f in recs[:30]:
                fi=f.get("id")
                filename=str(f.get("path") or f.get("filename") or f.get("name") or "")
                rawsize=f.get("size") or f.get("sizeBytes")
                z={"original_version_ID":ver,"original_file_ID":fi,
                   "original_file_name":filename,"declared_size":rawsize,
                   "download_status":"SOURCE_STOP"}
                if isinstance(fi,int) and filename and (rawsize is None or isinstance(rawsize,int) and 0<rawsize<=MAX_BYTES):
                    for attempt in (
                        f"{DRYAD_API}/files/{fi}/download",
                        f"https://datadryad.org/downloads/file_stream/{fi}",
                    ):
                        try:
                            bb=pull(attempt)
                            z["original_file_SHA256"]=hashlib.sha256(bb).hexdigest()
                            z["original_file_size_bytes"]=len(bb)
                            z["original_file_schema_only"]=safe_schema(filename,bb)
                            z["download_status"]="VERIFIED_RAW_BYTES"
                            if z["original_file_schema_only"]["kind"] in ("csv","tsv","xlsx","zip"):
                                out["source_authenticated_original_table_count"]+=1
                            break
                        except (ValueError,ImportError,zipfile.BadZipFile,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
                            z.setdefault("attempt_errors",[]).append({
                                "type":type(exc).__name__,
                                **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})})
                out["files"].append(z)
        except (ValueError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            out.setdefault("file_index_errors",[]).append({
                "version_ID":ver,"error_type":type(exc).__name__,
                **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})})
    out["metadata"]="VERIFIED_DRYAD_PUBLIC_API" if out.get("dataset_metadata_status")=="VERIFIED_DRYAD_PUBLIC_API" else "SOURCE_STOP"
    return out


def run():
    fig=get_figshare()
    dry=get_dryad()
    return {
        "schema":"alan_original_54_mesocosm_three_papers_figshare_dryad_source_gate_v1",
        "physical_programme":"Penryn_54_trough_2012_2016",
        "independent_light_experiment_programmes":1,
        "original_physical_light_mesocosms":54,
        "independent_mesocosms_per_light_trophic_cell":6,
        "crossed_light_arms":3,
        "crossed_trophic_arms":3,
        "original_publications":[
            {"doi":"10.1098/rstb.2014.0131","stage":"initial_plant_aphid_resource"},
            {"doi":"10.1111/1365-2664.13240","stage":"predator_herbivore_plant_resource"},
            {"doi":"10.1002/ece3.9157","stage":"long_term_plant_biomass_carbon_nitrogen"},
        ],
        "figshare_2022_plant_original_dataset":fig,
        "dryad_2018_herbivory_original_dataset":dry,
        "not_same_as_Lizard30_light_experiment":True,
        "source_54_plot_light_treatment_year_full_join_verified":False,
        "randomized_light_x_predator_effect_reanalyzed_from_raw":False,
        "cross_trophic_causal_mediation_estimated":False,
        "source_original_row_level_observations_exported":False,
        "source_original_plot_IDs_exported":False,
        "green_CI_does_not_imply_original_download":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=run()
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n")
    for key in ("figshare_2022_plant_original_dataset","dryad_2018_herbivory_original_dataset"):
        d=z[key]
        print("ALAN_54_PROGRAMME_PUBLIC_SOURCE",key,d["metadata"],"actual_tables",d["source_authenticated_original_table_count"])
        print("ALAN_54_PROGRAMME_SOURCE_FILES",key,d.get("files"))
    print("ALAN_54_MESOCOSM_3PUBLICATIONS_1EXPERIMENT_NO_NEW_META",z["source_54_plot_light_treatment_year_full_join_verified"])


if __name__=="__main__":
    main()
