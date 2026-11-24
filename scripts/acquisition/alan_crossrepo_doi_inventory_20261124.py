#!/usr/bin/env python3
"""Bounded, inspectable cross-repository DOI discovery for ALAN.

Scans the Git TREE PATHS of every public GitHub repository owned by USER,
then reads a deterministic priority subset of bibliography/source documents
from each repository's pinned default-branch tree.

This is a DISCOVERY audit, not 62 complete full-text reviews or a PRISMA
systematic database search. Unread candidate files and unindexed binary data
are explicitly counted. Never interpret absent DOI as an excluded study.

Output rows are paths, DOI identifiers and line numbers only; no copyrighted
source passages or original PDFs are redistributed.
"""
from __future__ import annotations

import argparse
import base64
import csv
import json
import os
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

USER="zuizui0223"
API="https://api.github.com"
DOI=re.compile(r"(?i)\b10\.\d{4,9}/[-._;()/:A-Z0-9]+")
READABLE={".md",".txt",".rst",".tex",".bib",".ris",".nbib",".enw",".csv",".tsv",".json"}
BIB={".bib",".ris",".nbib",".enw"}
SKIP_DIR={"/.git/","/.github/","/node_modules/","/vendor/","/__pycache__/","/venv/"}
MAX_SOURCE_BYTES=150_000
MAX_DOIS_PER_FILE=150
HEADERS={"Accept":"application/vnd.github+json",
         "User-Agent":"alan-scientific-crossrepo-doi-source-audit/2026-11-24",
         "X-GitHub-Api-Version":"2022-11-28"}
if os.environ.get("GITHUB_TOKEN"):
    HEADERS["Authorization"]="Bearer "+os.environ["GITHUB_TOKEN"]

def api(url,*,max_bytes=10_000_000,retry=2):
    if not url.startswith("https://api.github.com/"):
        raise ValueError("STOP: only official GitHub REST API permitted")
    for attempt in range(retry+1):
        try:
            with urlopen(Request(url,headers=HEADERS),timeout=25) as r:
                if int(r.headers.get("Content-Length") or "0")>max_bytes:
                    raise ValueError("STOP: GitHub payload too large")
                body=r.read(max_bytes+1)
            if len(body)>max_bytes:raise ValueError("STOP: downloaded body too large")
            return json.loads(body)
        except HTTPError as e:
            if e.code in (403,429,500,502,503,504) and attempt<retry:
                # Only bounded retry. Never claim completeness after failure.
                time.sleep(2*(attempt+1))
                continue
            raise
        except (TimeoutError,URLError):
            if attempt<retry:
                time.sleep(2*(attempt+1))
                continue
            raise

def list_repositories(owner):
    repos=[]
    for page in range(1,8):
        url=f"{API}/users/{quote(owner)}/repos?per_page=100&type=owner&page={page}"
        raw=api(url,max_bytes=8_000_000)
        if not isinstance(raw,list):raise ValueError("STOP: repo API response is not a list")
        for x in raw:
            if not x.get("name") or x.get("owner",{}).get("login","").lower()!=owner.lower():
                raise ValueError("STOP: source repo not owned by requested account")
            repos.append(x)
        if len(raw)<100:break
    else:
        raise ValueError("STOP: user has >700 repos; pagination bounds insufficient")
    ids=[r["full_name"].lower() for r in repos]
    if len(set(ids))!=len(ids):raise ValueError("STOP: duplicate repo IDs")
    return sorted(repos,key=lambda x:x["full_name"].lower())

def candidate_rank(path,byte_size):
    if byte_size is None or byte_size>MAX_SOURCE_BYTES or byte_size<1:
        return None
    lower=path.lower()
    if any(x in "/"+lower for x in SKIP_DIR):return None
    ext=Path(lower).suffix
    if ext not in READABLE:return None
    parts=lower.split("/")
    stem=Path(lower).stem
    score=1
    if ext in BIB:score+=190
    if re.search(r"(^|/)(references?|bibliograph(y|ies)|literature|citations?)([._/-]|$)",lower):score+=150
    if lower in ("readme.md","readme.txt"):score+=130
    if re.search(r"(alan|light.?pollut|nocturn|night.?light|pollinat|arthropod|insect|meta.analysis)",lower):score+=85
    if re.search(r"(source|primary|evidence|review|screen|experiment|stud(y|ies)|paper|publication)",lower):score+=50
    if "manuscript" in parts or "docs" in parts:score+=28
    if lower.endswith("references.bib"):score+=80
    if ext==".json" or ext==".csv":score-=6
    if "/results/" in "/"+lower or "/data/raw/" in "/"+lower:score-=40
    # Quality controls: do not accidentally choose tests and generated plots
    if lower.startswith("tests/") or lower.startswith("scripts/"):score-=50
    return score

def normalize_doi(raw):
    d=raw.strip().lower().rstrip(".,;:)")
    d=d.lstrip("([")
    # Preserve parentheses inside a DOI when balanced. End punctuation is
    # ambiguous; do not use guessed validity to call papers "included".
    return d

def search_document(repo,tree_path,blob_sha):
    uri=f"{API}/repos/{quote(repo, safe='/')}/git/blobs/{quote(blob_sha)}"
    blob=api(uri,max_bytes=650_000)
    if blob.get("sha")!=blob_sha:
        raise ValueError("STOP: current blob identity changed")
    encoding=blob.get("encoding")
    if encoding=="base64":
        data=base64.b64decode(blob.get("content",""),validate=False)
        if len(data)>MAX_SOURCE_BYTES:raise ValueError("STOP: blob too large after base64")
        txt=data.decode("utf-8",errors="replace")
    elif encoding=="utf-8":
        txt=blob.get("content","")
    else:raise ValueError("STOP: unsupported blob encoding")
    occurrences=[]
    seen=set()
    for line_number,line in enumerate(txt.splitlines(),start=1):
        if not ("10." in line or "doi.org" in line.lower()):continue
        for match in DOI.finditer(line):
            doi=normalize_doi(match.group())
            if len(doi)>175 or not doi:continue
            # Deduplicate the same DOI within a single source file, keep
            # first line for traceability, never count occurrences as studies.
            if doi not in seen:
                seen.add(doi)
                occurrences.append({"repo":repo,"path":tree_path,"git_blob_sha":blob_sha,
                                    "line":line_number,"doi":doi})
                if len(occurrences)>=MAX_DOIS_PER_FILE:
                    return occurrences,True
    return occurrences,False

def audit(owner,source_cap=8,boost=None):
    if boost is None:boost={"alan":40,"iwe":32,"mois":16}
    repos=list_repositories(owner)
    audit_rows=[]; citations=[]; source_failures=[]
    for index,repo in enumerate(repos,1):
        full=repo["full_name"];short=repo["name"].lower()
        branch=repo.get("default_branch") or "main"
        row={"repository":full,"default_branch":branch,"public_visibility":repo.get("visibility"),
             "archive_status":repo.get("archived"),"files_in_tree":"","tree_sha":"",
             "tree_complete":False,"candidate_text_files":0,"candidate_bibliography_files":0,
             "candidate_bytes_over_limit":0,"target_source_files":0,"read_source_files":0,
             "unread_priority_files":0,"document_content_scan_complete":False,
             "new_doI_literals_found":0,"source_read_errors":0,"truncated_file_doi_occurrences":0,
             "audit_status":"STARTED","error":""}
        try:
            tree=api(f"{API}/repos/{quote(full,safe='/')}/git/trees/{quote(branch,safe='')}?recursive=1",max_bytes=20_000_000)
            files=[x for x in tree.get("tree",[]) if x.get("type")=="blob"]
            row["files_in_tree"]=len(files)
            row["tree_sha"]=tree.get("sha","")
            row["tree_complete"]=not bool(tree.get("truncated"))
            row["candidate_bytes_over_limit"]=sum(1 for x in files if Path(x["path"].lower()).suffix in READABLE and (x.get("size") or 0)>MAX_SOURCE_BYTES)
            ranked=[]
            for obj in files:
                p=obj["path"]
                value=candidate_rank(p,obj.get("size"))
                if value is not None:ranked.append((value,p,obj))
            ranked.sort(key=lambda x:(-x[0],x[1]))
            row["candidate_text_files"]=len(ranked)
            row["candidate_bibliography_files"]=sum(1 for _,p,_ in ranked if Path(p.lower()).suffix in BIB)
            cap=boost.get(short,source_cap)
            selected=ranked[:cap]
            row["target_source_files"]=len(selected)
            row["unread_priority_files"]=max(0,len(ranked)-len(selected))
            for _,path,obj in selected:
                try:
                    found,more=search_document(full,path,obj["sha"])
                    citations.extend(found)
                    row["read_source_files"]+=1
                    row["truncated_file_doi_occurrences"]+=int(more)
                except (Exception) as exc:
                    row["source_read_errors"]+=1
                    source_failures.append({"repository":full,"path":path,"kind":type(exc).__name__,
                                            "error":str(exc)[:180]})
            row["new_doI_literals_found"]=len(set(c["doi"] for c in citations if c["repo"]==full))
            row["document_content_scan_complete"]=bool(row["tree_complete"] and row["unread_priority_files"]==0 and row["source_read_errors"]==0 and row["candidate_bytes_over_limit"]==0 and row["truncated_file_doi_occurrences"]==0)
            row["audit_status"]="SAMPLED_DOI_SOURCES" if row["read_source_files"] else "PATH_ONLY"
        except Exception as exc:
            row["audit_status"]="TREE_ACCESS_STOP"
            row["error"]=f"{type(exc).__name__}: {str(exc)[:220]}"
        audit_rows.append(row)
        print(f"CROSSREPO_PROGRESS {index}/{len(repos)} {full} tree={row['tree_complete']} source={row['read_source_files']}/{row['candidate_text_files']} doi={row['new_doI_literals_found']} status={row['audit_status']}",flush=True)

    grouped={}
    for x in citations:
        z=grouped.setdefault(x["doi"],{"doi":x["doi"],"repositories":set(),"sources":0})
        z["repositories"].add(x["repo"]);z["sources"]+=1
    doi_summaries=[{"doi":x["doi"],"repo_count":len(x["repositories"]),
                    "source_files":x["sources"],"repositories":";".join(sorted(x["repositories"]))}
                   for x in grouped.values()]
    doi_summaries.sort(key=lambda x:(-x["repo_count"],x["doi"]))
    return audit_rows,citations,doi_summaries,source_failures

def save_csv(file,rows,fields):
    with file.open("w",encoding="utf-8",newline="") as out:
        w=csv.DictWriter(out,fieldnames=fields)
        w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--owner",default=USER)
    ap.add_argument("--output-dir",type=Path,default=Path("results/alan_crossrepo_doi_inventory_20261124"))
    ap.add_argument("--source-cap-per-repo",type=int,default=8)
    args=ap.parse_args()
    if args.owner!=USER:raise ValueError("STOP: audit owner differs from predeclared source account")
    if not 1<=args.source_cap_per_repo<=50:
        raise ValueError("STOP: invalid bounded file budget")
    coverage,citations,dois,failures=audit(args.owner,source_cap=args.source_cap_per_repo)
    p=args.output_dir;p.mkdir(parents=True,exist_ok=True)
    save_csv(p/"repository_coverage.csv",coverage,list(coverage[0]))
    save_csv(p/"doi_occurrence_sources.csv",citations,["repo","path","git_blob_sha","line","doi"])
    save_csv(p/"doi_crossrepo_duplicates.csv",dois,["doi","repo_count","source_files","repositories"])
    save_csv(p/"source_access_failures.csv",failures,["repository","path","kind","error"])
    n=len(coverage);sum_col=lambda col:sum(int(x[col]) for x in coverage)
    summary={
       "date":"2026-11-24","scope":"Default branches of all owned public GitHub repositories; path-level exhaustive where recursive tree nontruncated",
       "repositories_reported":n,
       "repositories_with_full_tree":sum(bool(x["tree_complete"]) for x in coverage),
       "repositories_path_audit_failed":[x["repository"] for x in coverage if not x["tree_complete"]],
       "tracked_file_entries":sum_col("files_in_tree"),
       "candidate_text_files":sum_col("candidate_text_files"),
       "priority_files_actually_read":sum_col("read_source_files"),
       "priority_candidate_unread":sum_col("unread_priority_files"),
       "bibliography_files_path_identified":sum_col("candidate_bibliography_files"),
       "DOI_source_occurrences_distinct_by_file":len(citations),
       "distinct_literal_DOIs_in_read_sources":len(dois),
       "DOIs_appearing_in_multiple_repositories":sum(x["repo_count"]>1 for x in dois),
       "coverage_is_full_text_systematic_review":False,
       "all_default_branch_text_files_scanned":all(x["document_content_scan_complete"] for x in coverage),
       "all_references_full_text_adjudicated":False,
       "unmerged_branches_PRs_and_Git_history_covered":False,
       "DOI_matches_primary_ALAN_study_count":None,
       "final_PRISMA_eligible_study_count":None,
       "research_caveat":"DOI appearances include background citations, reviews, code comments and reused references; they are not independent source experiments. Do not add to meta k.",
    }
    (p/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("CROSSREPO_FINAL",json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
