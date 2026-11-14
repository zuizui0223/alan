#!/usr/bin/env python3
"""Discover corrected 2025 Dryad moth light-activity experimental source.

Original Briolat et al. (2026 Proc B) DOI 10.1098/rspb.2025.2704.
Dataset DOI 10.5061/dryad.jq2bvq8ps version updated Dec 11 2025.
This first phase retrieves only full public *dataset metadata and files
manifests*, not video frame rows. No new activity effect is inferred here.
"""
from __future__ import annotations
import argparse,json,urllib.parse,urllib.request,urllib.error
from pathlib import Path
BASE="https://datadryad.org/api/v2"
DOI="doi:10.5061/dryad.jq2bvq8ps"
ENC=urllib.parse.quote(DOI,safe="")
def req(url,timeout=30):
    request=urllib.request.Request(url,headers={
        "Accept":"application/json","User-Agent":"ALAN-scientific-source-audit/1.0"})
    with urllib.request.urlopen(request,timeout=timeout) as r:
        return json.load(r)
def hrefs(obj):
    if isinstance(obj,dict):
        for k,v in obj.items():
            if k=="href" and isinstance(v,str):
                yield v
            else:
                yield from hrefs(v)
    elif isinstance(obj,list):
        for x in obj:yield from hrefs(x)
def objects(obj):
    if isinstance(obj,dict):
        yield obj
        for v in obj.values():yield from objects(v)
    elif isinstance(obj,list):
        for x in obj:yield from objects(x)
def safe_file(o):
    return {k:o.get(k) for k in ("id","path","name","size","digest","digestType",
           "mimeType","status","fileType","description") if k in o}
def main(out):
    out=Path(out)
    if out.exists():raise ValueError("Refuse overwrite")
    url=BASE+"/datasets/"+ENC
    original=req(url)
    versions=req(url+"/versions")
    v=[]
    for obj in objects(versions):
        if (("versionNumber" in obj or "version" in obj)
                and ("id" in obj or "versionNumber" in obj)):
            v.append({k:obj[k] for k in ("id","versionNumber","version","versionStatus",
                "publicationDate","lastModified") if k in obj})
    paths=[x for x in hrefs(versions) if "/versions/" in x or "/files" in x]
    version_debug=[{"keys":sorted(o.keys()),"versionNumber":o.get("versionNumber"),
                    "links":o.get("_links")} for o in objects(versions)
                   if isinstance(o,dict) and "versionNumber" in o][:12]
    # Follow only known origin; do not allow arbitrary redirected external APIs.
    for x in paths:
        if not x.startswith("https://"):
            continue
        if not x.startswith("https://datadryad.org/"):
            continue
        # Version file lists use /api/v2/versions/<id>/files.
    # Dryad public V2 returns relative links and numeric version IDs only in
    # each version self URL; a version NUMBER (4) is not a version ID (413461).
    source_versions=[
       o for o in objects(versions) if isinstance(o,dict) and
       isinstance(o.get("versionNumber"),int) and
       isinstance(o.get("_links"),dict)]
    if not source_versions:raise ValueError("No published dataset versions")
    selected=max(source_versions,key=lambda v:v["versionNumber"])
    if selected["versionNumber"] != 4:
        raise ValueError("Corrected 2025-12-11 Dryad version changed")
    href=selected["_links"].get("stash:files",{}).get("href")
    if not isinstance(href,str):raise ValueError("No corrected source files link")
    candidate=urllib.parse.urljoin("https://datadryad.org",href)
    if not candidate.startswith(BASE+"/versions/"):
        raise ValueError("Unexpected untrusted Dryad corrected source origin")
    candidates=[candidate]
    latest_public_version=selected["versionNumber"]
    errors=[]
    files=None
    for candidate in candidates:
        try:
            z=req(candidate)
            potential=[obj for obj in objects(z) if
                    ("path" in obj or "name" in obj) and
                    ("size" in obj or "digest" in obj)]
            if potential:
                files=[safe_file(x) for x in potential]
                break
        except (urllib.error.HTTPError,urllib.error.URLError,ValueError) as e:
            errors.append({"target":candidate,"error":type(e).__name__,
                           "http_status":getattr(e,"code",None)})
    result={"schema":"ALAN_Dryad_original_moth_23species_activity_public_file_manifest_v1",
        "status":"ORIGINAL_PUBLIC_DATASET_MANIFEST_NO_ACTIVITY_RESULT",
        "DOI":DOI,
        "source_publication_DOI":"10.1098/rspb.2025.2704",
        "api_dataset":url,
        "dataset_keys":sorted(original),
        "source_dataset_current_version":original.get("mostRecentVersion"),
        "raw_versions_metadata_schema_keys":sorted(versions),
        "raw_versions_indexed":v,
        "listed_original_files":files,
        "corrected_selected_version_number":latest_public_version,
        "file_get_errors":errors,
        "data_not_downloaded_this_phase":True,
        "corrected_2025_12_11_data_preferred_over_2025_12_03":True,
        "main_843_moth_23species_three_families_numbers_from_published_abstract_not_recomputed":True,
        "no_video_frame_pseudoreplication_or_pollen_delivery_claim":True,
        "url_candidates":candidates,
        "version_public_metadata_debug":version_debug,
        "version_link_href_samples":paths[:15]}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_DRYAD_MOTH_23SP_MANIFEST",json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args(); main(a.out)
