#!/usr/bin/env python3
"""Verify Xue et al. 2026 original LIVE macrophyte ALAN x snail source.

Paper: 10.1002/ecy.70458
Original author-deposited study: 10.57760/sciencedb.35715
Publisher's data-availability statement points to Science DataBank; the
repository's official 2026 public OpenAPI declares GET /json?doi=...
(https://www.scidb.cn/open-api/v2/api-docs).
Do not substitute figure-digitization, guesses or published percent means
for the actual per mesocosm/zone source rows.

Original designs (article):
 - monoculture: 6 independently exposed zones (3 lit / 3 dark),
   4 plant taxa, 80 cylinders with 77 survivors; cylinders nested in zone.
 - feeding choice / no choice: 4-day original bucket assay; each source
   snail x plant x light combination n=5 as described by authors, but
   source light fixture allocation must be verified before inference.
 - mixed community: 2 x 3 design, 60 mesocosms, 57 final with biomass;
   physical light zone / shared fixture units are not identified in abstract.
Stop experimental light pseudo-replication when no zone IDs are provided.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error,urllib.parse,zipfile
from pathlib import Path

SCHEMA="alan_Xue2026_ScienceDB_35715_original_source_and_experimental_unit_gate_v1"
DOI="10.57760/sciencedb.35715"
SOURCE_ID="b9b02a9ec07c456b8584dd5fc732b7da"
PAPER="10.1002/ecy.70458"
MAX_META=4<<20
MAX_FILE=12<<20
URLS=[
 f"https://www.scidb.cn/api/sdb-openapi-service/json?doi={urllib.parse.quote(DOI,safe='')}",
 f"https://www.scidb.cn/open-api/v2/json?doi={urllib.parse.quote(DOI,safe='')}",
 f"https://www.scidb.cn/en/detail?dataSetId={SOURCE_ID}",
]
def request(url, maxbytes=MAX_META):
    req=urllib.request.Request(url, headers={
      "User-Agent":"alan-macro-ecology-source-reproduction/2026",
      "Accept":"application/json,text/html,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=35) as res:
        raw=res.read(maxbytes+1)
        if len(raw)>maxbytes:raise RuntimeError("SOURCE_TOO_LARGE_STOP")
        return raw,res.geturl(),res.headers.get("Content-Type","")

def scan_meta(node,depth=0):
    """Only browse metadata, never invent download IDs or raw measured columns."""
    if depth>12:return []
    files=[]
    if isinstance(node,dict):
        keys={str(k).lower():k for k in node}
        name=next((node[keys[k]] for k in ("filename","filenameori","file_name","name","key")
                   if k in keys and isinstance(node[keys[k]],str)),None)
        if name and (re.search(r"\.(zip|csv|xlsx|xls|tsv|rds|rda|txt)$",name,re.I)):
            links={str(k):v for k,v in node.items()
                   if isinstance(v,str) and
                   ("url" in str(k).lower() or "download" in str(k).lower()
                    or "path" in str(k).lower() or "link" in str(k).lower())}
            files.append({"reported_name":name[:220],
                  "metadata_size":next((node[keys[k]] for k in ("size","filesize","bytes")
                                        if k in keys),None),
                  "metadata_direct_links":links})
        for value in node.values():files+=scan_meta(value,depth+1)
    elif isinstance(node,list):
        for value in node[:5000]:files+=scan_meta(value,depth+1)
    return files

def source_schema(buf,name):
    original_sha=hashlib.sha256(buf).hexdigest()
    if name.lower().endswith((".csv",".tsv")):
        decoded=buf.decode("utf-8-sig",errors="replace")
        for sep in (",",";","\t"):
            rr=list(csv.DictReader(io.StringIO(decoded),delimiter=sep))
            if rr and len(rr[0])>=3:
                head=list(rr[0])
                return {"type":"csv_source","SHA256":original_sha,
                        "rows":len(rr),"columns":head[:120],
                        "source_factor_levels":{
                         col:sorted({r[col].strip() for r in rr})[:35]
                         for col in head
                         if any(x in col.lower() for x in ("light","alan","zone","species","snail","mesocosm"))
                         and len({r[col] for r in rr})<=35}}
        return {"type":"unrecognized_original_CSV","SHA256":original_sha}
    if name.lower().endswith((".xlsx",".zip")) and zipfile.is_zipfile(io.BytesIO(buf)):
        with zipfile.ZipFile(io.BytesIO(buf)) as z:
            if z.testzip() is not None:raise ValueError("original ZIP CRC failed")
            names=[q.filename for q in z.infolist()]
        return {"type":"verified_ZIP_or_XLSX","SHA256":original_sha,
                "members":len(names),"original_member_names":names[:60]}
    return {"type":"unsupported_source_media","SHA256":original_sha}

def execute():
    result={"schema":SCHEMA,"paper_DOI":PAPER,"original_dataset_DOI":DOI,
      "original_ScienceDB_ID_from_DOI":SOURCE_ID,
      "official_openapi_documentation":"https://www.scidb.cn/open-api/v2/api-docs",
      "source_attempts":[],
      "source_original_data_not_fabricated":True,
      "published_experiment_1_zone_count":6,
      "published_experiment_1_zones_per_light":3,
      "published_experiment_1_pot_total":80,
      "published_experiment_1_pots_after_mortality":77,
      "published_experiment_3_mesocosms":60,
      "published_experiment_3_biomass_complete_mesocosms":57,
      "source_light_zone_identifiers_verified":False,
      "source_species_level_light_by_herbivory_effect_fitted":False,
      "wet_affinity_trait_moderator_fitted":False,
      "source_original_experiment_level_variances_obtained":False,
      "source_files":[],"status":"ORIGINAL_SOURCE_ACCESS_STOP"}
    for u in URLS:
        record={"url":u}
        try:
            blob,final,ctype=request(u)
            record.update(status="HTTP_OK",bytes=len(blob),
                          final_url=final,mime=ctype,
                          source_bytes_SHA256=hashlib.sha256(blob).hexdigest())
            if "json" in (ctype or "").lower() or blob[:1] in (b"{",b"["):
                m=json.loads(blob)
                rep=json.dumps(m,ensure_ascii=False)
                record["data_record_DOI_literal_present"]=DOI.lower() in rep.lower()
                result["source_files"]+=scan_meta(m)
                record["metadata_top_level_keys"]=(list(m)[:50] if isinstance(m,dict) else [])
                result["status"]="SOURCE_METADATA_VERIFIED" if (
                  record["data_record_DOI_literal_present"]) else "SOURCE_METADATA_UNVERIFIED_STOP"
                # Stop at first genuinely matching metadata. Alternate page
                # URL is only a navigational status check, not source provenance.
                if record["data_record_DOI_literal_present"]:
                    result["source_metadata_SHA256"]=record["source_bytes_SHA256"]
                    break
            else:
                record["html_not_original_data"]=True
                # No page scraping guessed file URLs.
        except Exception as e:
            record.update(status="SOURCE_ACCESS_STOP",error_type=type(e).__name__,
                          http_code=getattr(e,"code",None),
                          error_excerpt=str(e)[:150])
        finally:
            result["source_attempts"].append(record)
    # Scan only links explicitly in genuine metadata; never guess numeric
    # download routes. Require source link from official ScienceDB domain.
    result["source_files"]=(list({f["reported_name"]:f for f in result["source_files"]}.values()))
    if result["status"]=="SOURCE_METADATA_VERIFIED":
        for item in result["source_files"]:
            metaurl=next((value for value in item["metadata_direct_links"].values()
                           if isinstance(value,str) and value.startswith("https://")
                           and (urllib.parse.urlparse(value).hostname or "").endswith("scidb.cn")),None)
            if not metaurl:
                item["status"]="NO_ORIGINAL_DOWNLOAD_LINK_METADATA_STOP"
                continue
            try:
                raw,final,ctype=request(metaurl,MAX_FILE)
                item["source_byte_verification"]=source_schema(raw,item["reported_name"])
                item["source_download_final_url"]=final
                item["status"]="ORIGINAL_FILE_SOURCE_VERIFIED"
            except Exception as e:
                item["status"]="ORIGINAL_DOWNLOAD_STOP"
                item["download_error_type"]=type(e).__name__
                item["download_http_code"]=getattr(e,"code",None)
        if any(x["status"]=="ORIGINAL_FILE_SOURCE_VERIFIED" for x in result["source_files"]):
            result["status"]="ORIGINAL_SOURCE_FILE_GENUINE_AVAILABLE_NEEDS_UNIT_QA"
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    z=execute()
    target=Path(a.out);target.parent.mkdir(exist_ok=True,parents=True)
    target.write_text(json.dumps(z,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("XUE2026_LIVE_PLANT_SCIENCEDB_SOURCE_GATE "+json.dumps(z,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
