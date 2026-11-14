#!/usr/bin/env python3
"""Audit the real Macgregor 2019 archived script's needed primary CSV files.

This only reads the author's public analysis code and repository file manifests,
not row-level reproductive outcomes. Assert exact archived blob identity and
record which source files are needed for a real statistical refit.
"""
import argparse
import hashlib
import json
import re
import urllib.error
import urllib.request
from pathlib import Path

AUTHOR="CallumJMacgregor/BioassayPaper"
SCRIPT="FinalScripts.R"
EXPECTED_SCRIPT_GIT_BLOB_SHA1="04136ea82d6c1bbef22f4fd60c41de9907b00550"
API_TREE=f"https://api.github.com/repos/{AUTHOR}/git/trees/master?recursive=1"
RAW=f"https://raw.githubusercontent.com/{AUTHOR}/master/{SCRIPT}"
ZENODO=f"https://zenodo.org/api/records/1743917"
UA={"User-Agent":"ALAN-Macgregor2019-source-data-citable-audit/1.0",
    "Accept":"application/json"}
NEEDED={"SeedSetBinom.csv","SeedSetYN.csv","SeedCount.csv"}

def require(x,msg):
    if not x: raise RuntimeError("ALAN_MACGREGOR_REAL2019: "+msg)

def grab(url,limit=2_000_000):
    req=urllib.request.Request(url,headers=UA)
    with urllib.request.urlopen(req,timeout=50) as r:
        b=r.read(limit+1)
    require(len(b)<=limit,"File oversize or truncated")
    return b

def main(out):
    out=Path(out); require(not out.exists(),"Do not overwrite result")
    code=grab(RAW)
    blob=hashlib.sha1(b"blob "+str(len(code)).encode()+b"\0"+code).hexdigest()
    require(blob==EXPECTED_SCRIPT_GIT_BLOB_SHA1,
        "Published code blob identity changed: "+blob)
    s=code.decode("utf-8-sig")
    # Original Windows path written in R string as Data\\filename.csv.
    cited=sorted(set(re.findall(r"Seed(?:Set(?:Binom|YN)|Count)\.csv",s,re.I)))
    require(set(cited)==NEEDED,"Author code input references changed: "+str(cited))
    tree=json.loads(grab(API_TREE))
    paths={n.get("path") for n in tree.get("tree",[])}
    require(SCRIPT in paths,"Author original R script not in repository")
    # Zenodo archive may be source code only; verify archive members explicitly
    # via source self link if available rather than guessing from its small size.
    metadata=json.loads(grab(ZENODO))
    require(int(metadata.get("id",-1))==1743917,"Zenodo ID mismatch")
    files=metadata.get("files",[])
    manifest=[{"name":f.get("key") or f.get("filename"),
        "size":f.get("size"),"checksum":f.get("checksum")} for f in files]
    # No original data file in the GitHub tree and no raw csv registered
    # as its own Zenodo artifact does NOT prove raw data not inside zip.
    archive_member_names=[]
    archive_result="NOT_DOWNLOADED"
    potential=[f for f in files if str(f.get("key","")).lower().endswith(".zip")]
    if len(potential)==1:
        target=potential[0]
        link=target.get("links",{}).get("self") or target.get("links",{}).get("download")
        if link:
            import io,zipfile
            try:
                raw=grab(link,limit=2_000_000)
                if zipfile.is_zipfile(io.BytesIO(raw)):
                    with zipfile.ZipFile(io.BytesIO(raw)) as zipf:
                        archive_member_names=zipf.namelist()
                    archive_result="AUTHOR_ZENODO_ZIP_MEMBERS_READ"
                else:
                    archive_result="ZENODO_RETURNED_NON_ZIP_NO_RAW_CLAIM"
            except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError):
                archive_result="AUTHOR_ZENODO_ZIP_ACCESS_BLOCKED"
        else: archive_result="METADATA_NO_DOWNLOAD_LINK"
    needed_found_repo=sorted(p for p in paths if p.rsplit("/",1)[-1] in NEEDED)
    needed_found_zip=sorted(p for p in archive_member_names if p.rsplit("/",1)[-1] in NEEDED)
    result={"schema":"ALAN_Macgregor2019_real_archived_analysis_code_source_input_availability_v1",
        "author_git_repo":AUTHOR,"author_source_script":SCRIPT,
        "exact_author_script_sha1_verified":blob,
        "required_original_primary_csv_filenames":cited,
        "original_required_csv_present_in_public_GitHub_tree":needed_found_repo,
        "Zenodo_DOI":"10.5281/zenodo.1743917",
        "zenodo_record_files":manifest,
        "original_zip_member_listing_status":archive_result,
        "original_zip_members":archive_member_names,
        "required_original_CSV_names_found_in_Zenodo_zip":needed_found_zip,
        "2019_full_numerical_GLM_reproduced":False,
        "negative_source_claim_qualified":"No original primary CSV was available in scanned GitHub tree; Zenodo ZIP status separately captured",
        "rule":"An archived and SHA verified analysis script is NOT an observed numeric dataset; never synthesize plant observations or a published effect-size SE."}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_MACGREGOR2019_CODE_SOURCE_DEPENDENCY",json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args().out)
