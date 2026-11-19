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


def pdf_schema(source):
    """Read original archival PDF text/table geometry; never publish raw cells."""
    from pypdf import PdfReader
    if not source.startswith(b"%PDF-"):
        raise ValueError("Not a real published PDF source")
    reader=PdfReader(io.BytesIO(source))
    if len(reader.pages)>80:
        raise ValueError("Unusually long original data PDF; inspect separately")
    pages=[]
    for p in reader.pages:
        text=p.extract_text() or ""
        lines=text.splitlines()
        # Numeric-heavy rows provide evidence of tabulated individual data,
        # not actual values. Do not export source plot or treatment rows.
        numeric_lines=sum(bool(re.findall(r"-?\d+(?:\.\d+)?",line))
                          and len(re.findall(r"-?\d+(?:\.\d+)?",line))>=4
                          for line in lines)
        pages.append({
            "extractable_characters":len(text),
            "lines_with_four_or_more_numeric_fields":numeric_lines,
            "page_line_count":len(lines),
            "contains_mesocosm_keyword":bool(re.search(r"mesocosm|trough",text,re.I)),
            "contains_biomass_keyword":bool(re.search(r"biomass|dry\s*mass",text,re.I)),
            "contains_treatment_keyword":bool(re.search(r"treatment|light|predator|herbivor",text,re.I)),
        })
    return {"kind":"pdf","page_count":len(reader.pages),
            "page_text_shapes_only":pages,
            "no_original_pdf_measurement_cells_exported":True,
            "cannot_claim_plant_original_numeric_table_reanalysis_from_PDF_shape_only":True}


def safe_schema(name,b):
    suffix=Path(name).suffix.lower()
    if suffix==".csv":
        return {"kind":"csv","sheets":[csv_schema(b)]}
    if suffix==".tsv":
        return {"kind":"tsv","sheets":[csv_schema(b,delimiter="\t")]}
    if suffix==".xlsx":
        return {"kind":"xlsx","sheets":xlsx_schema(b)}
    if suffix==".pdf":
        return pdf_schema(b)
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


def dryad_official_href(href, *, allow_files=False):
    """Follow the HAL API link, not the absent integer 'id' metadata field."""
    if not isinstance(href,str) or not href:
        raise ValueError("Original Dryad HAL link missing")
    url=urllib.parse.urljoin("https://datadryad.org",href)
    parts=urllib.parse.urlsplit(url)
    if parts.scheme!="https" or parts.hostname!="datadryad.org" or parts.port not in (None,443):
        raise ValueError("Refuse nonofficial Dryad HAL source")
    path=parts.path
    if not path.startswith("/api/v2/"):
        raise ValueError("Refuse guessed or non-API Dryad link")
    if not allow_files and not re.fullmatch(r"/api/v2/versions/[0-9]+/files",path):
        raise ValueError("Link is not a genuine original Dryad version file list")
    if allow_files and not (re.fullmatch(r"/api/v2/files/[0-9]+/download",path) or
                            re.fullmatch(r"/api/v2/versions/[0-9]+/files",path)):
        raise ValueError("Link is not a genuine original Dryad file download")
    return url


def dryad_hal_version_file_lists(dataset,versions):
    """Return the true official /versions/<id>/files URLs from HAL links.

    Older Dryad JSON has _embedded.stash:versions with _links.stash:files,
    but does not expose an integer 'id'. Extracting rec['id'] silently
    yielded zero file lists, which falsely looked like inaccessible raw data.
    """
    candidates=[]
    records=(versions.get("_embedded") or {}).get("stash:versions",[])
    if not isinstance(records,list):
        raise ValueError("Unrecognized Dryad versions HAL collection")
    for rec in records:
        if not isinstance(rec,dict):
            continue
        links=rec.get("_links") or {}
        file_link=(links.get("stash:files") or {}).get("href")
        if file_link:
            candidates.append(dryad_official_href(file_link))
        elif isinstance(rec.get("id"),int):
            candidates.append(dryad_official_href(f"/api/v2/versions/{rec['id']}/files"))
    cur=(dataset.get("_links") or {}).get("stash:version") or {}
    verhref=cur.get("href","")
    match=re.search(r"^/api/v2/versions/([0-9]+)$",urllib.parse.urlsplit(verhref).path)
    if match:
        candidates.append(dryad_official_href(f"/api/v2/versions/{match.group(1)}/files"))
    return list(dict.fromkeys(candidates))


def get_dryad():
    out={"dataset_DOI":DRYAD_DOI,"metadata":"SOURCE_STOP",
         "files":[],"source_authenticated_original_table_count":0,
         "original_light_meso_ids_jointly_matched_other_source":False}
    doi=urllib.parse.quote("doi:"+DRYAD_DOI,safe="")
    source_dataset={}
    source_versions={}
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
                source_dataset=data
                out["dataset_title"]=data.get("title","")
                out["dataset_version_metadata_keys"]=sorted(data.keys())
            if label=="versions":
                source_versions=data
                out["versions_metadata_keys"]=sorted(data.keys()) if isinstance(data,dict) else []
                embedded=data.get("_embedded",{}) if isinstance(data,dict) else {}
                ids=[]
                for v in embedded.values():
                    if isinstance(v,list):
                        for rec in v:
                            if isinstance(rec,dict) and isinstance(rec.get("id"),int):
                                ids.append(rec["id"])
                out["resolved_version_IDs"]=ids[:5]
                out["version_first_layer"]={
                  "top_level_keys":list(data.keys())[:30],
                  "embedded_collection_keys":list(embedded.keys()),
                  "embedded_shapes":{
                    str(k):{"value_type":type(v).__name__,
                            "n_entries":len(v) if isinstance(v,list) else None,
                            "first_record_keys":list(v[0].keys()) if isinstance(v,list) and v and isinstance(v[0],dict) else []}
                    for k,v in embedded.items()
                  },
                }
                # Dryad older deposits sometimes expose their current
                # version through a version "href" but use noninteger ids.
                # Do NOT guess file identifiers from DOI suffix or filenames.

        except (ValueError,KeyError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            out[label+"_metadata_status"]="SOURCE_STOP"
            out[label+"_error_type"]=type(exc).__name__
            if isinstance(exc,urllib.error.HTTPError):
                out[label+"_http_status"]=exc.code
    # The source API returned versions without a numeric ID. The *actual*
    # file-list path is in each HAL record's _links.stash:files.href.
    resolved_urls=dryad_hal_version_file_lists(source_dataset,source_versions)
    out["verified_original_version_file_list_urls"]=resolved_urls
    out["resolved_via_original_HAL_links"]=bool(resolved_urls)
    for file_list_url in resolved_urls[:4]:
        ver_match=re.search(r"/versions/([0-9]+)/files",file_list_url)
        ver=int(ver_match.group(1)) if ver_match else None

        try:
            blob=pull(file_list_url,limit=2_000_000)
            doc=json.loads(blob)
            out.setdefault("verified_versions_file_index",[]).append({
                "version_ID":ver,
                "original_HAL_file_list_URL":file_list_url,
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
                if not isinstance(fi,int):
                    links=f.get("_links") or {}
                    dl=(links.get("stash:download") or {}).get("href","")
                    match=re.search(r"/files/([0-9]+)/download",urllib.parse.urlsplit(dl).path)
                    if match: fi=int(match.group(1))
                filename=str(f.get("path") or f.get("filename") or f.get("name") or "")
                rawsize=f.get("size") or f.get("sizeBytes")
                z={"original_version_ID":ver,"original_file_ID":fi,
                   "original_file_name":filename,"declared_size":rawsize,
                   "download_status":"SOURCE_STOP"}
                dl=(f.get("_links") or {}).get("stash:download") or {}
                official_download=[]
                if isinstance(dl.get("href"),str):
                    try:
                        official_download.append(dryad_official_href(dl["href"],allow_files=True))
                    except ValueError:
                        z["source_download_link_status"]="HAL_SOURCE_UNEXPECTED"
                if isinstance(fi,int):
                    official_download.append(f"{DRYAD_API}/files/{fi}/download")
                    # These two public browser links are documented by
                    # Dryad's download UI, distinct from the API's
                    # sometimes authenticated file endpoint.
                    official_download.extend((
                        f"https://datadryad.org/downloads/file_stream/{fi}",
                        f"https://datadryad.org/stash/downloads/file_stream/{fi}",
                    ))
                if filename and (rawsize is None or isinstance(rawsize,int) and 0<rawsize<=MAX_BYTES):
                    for attempt in list(dict.fromkeys(official_download)):
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
    # Additional official Dryad archive endpoint documented in API:
    # GET /datasets/{doi}/download. We do not infer raw data availability
    # from metadata alone; only verified ZIP member schemas count.
    if out["source_authenticated_original_table_count"]==0:
        verfileid=re.search(r"/versions/([0-9]+)/files",resolved_urls[0]) if resolved_urls else None
        routes=[("original_dataset_archive",f"{DRYAD_API}/datasets/{doi}/download")]
        if verfileid:
            routes.append(("original_version_archive",
                           f"{DRYAD_API}/versions/{verfileid.group(1)}/download"))
        for label,url in routes:
            try:
                b=pull(url,limit=MAX_BYTES)
                if not b.startswith(b"PK\x03\x04"):
                    raise ValueError("Official Dryad archive is not original ZIP bytes")
                schema=safe_schema("original_dataset.zip",b)
                if not schema["members"]:
                    raise ValueError("Original Dryad archive has no tabular source files")
                out["official_archive_status"]="ORIGINAL_SOURCE_ARCHIVE_ZIP_VERIFIED"
                out["official_archive_SHA256"]=hashlib.sha256(b).hexdigest()
                out["official_archive_size_bytes"]=len(b)
                out["official_archive_schema_only"]=schema
                out["source_authenticated_original_table_count"]=len(schema["members"])
                out["official_archive_source_URL"]=url
                break
            except (ValueError,zipfile.BadZipFile,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
                out.setdefault("official_archive_attempts",[]).append({
                    "source_route":label,
                    "status":"SOURCE_STOP","error_type":type(exc).__name__,
                    **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})
                })
    out["retrieval_status"]=(
        "RAW_ORIGINAL_TABLES_VERIFIED" if out["source_authenticated_original_table_count"]>0
        else "FILE_METADATA_ONLY_OR_SOURCE_STOP"
    )
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
        if key.startswith("dryad") and d.get("official_archive_status")=="ORIGINAL_SOURCE_ARCHIVE_ZIP_VERIFIED":
            archive=d["official_archive_schema_only"]
            print("ALAN_DRYAD_ORIGINAL_DATASET_ARCHIVE_SHA",
                  d["official_archive_SHA256"])
            print("ALAN_DRYAD_ORIGINAL_DATASET_CSV_SCHEMAS",
                  [{"member":entry.get("archive_member"),
                    "sha":entry.get("member_SHA256"),
                    "sheets":entry.get("sheets",[])}
                   for entry in archive.get("members",[])])

        if key.startswith("dryad"):
            print("ALAN_54_DRYAD_SOURCE_METADATA_DIAG",{k:v for k,v in d.items() if k.endswith("_status") or k.endswith("_keys") or k in ("resolved_version_IDs","version_first_layer","verified_versions_file_index","file_index_errors","dataset_title")})
    print("ALAN_54_MESOCOSM_3PUBLICATIONS_1EXPERIMENT_NO_NEW_META",z["source_54_plot_light_treatment_year_full_join_verified"])


if __name__=="__main__":
    main()
