#!/usr/bin/env python3
"""Public 2017 moth feeding ALAN sensitivity; pathway gate, NOT plant fertility.

van Langevelde et al. 2017, Biol Lett doi:10.1098/rsbl.2016.0874.
Public CC0 Dryad source doi:10.5061/dryad.8v212, Data_Moth_Foraging.txt.

Original observations every ~6 minutes are NOT independent moths.
Summarise individual per-source-unique.id feeding fractions, then light×
species block summaries. Never use as pollen transfer, seed success, or
same-population plant nocturnal dependence V.

Runs official Dryad v2 anonymous API first; last-resort public legacy
file_stream is attempted without bypassing access controls. Fail closed
with detailed HTTP receipts if file bytes cannot be obtained. Do not
commit original individual data; aggregate report only.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import io
import json
import math
import pathlib
import statistics
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime,timezone

DOI="10.5061/dryad.8v212"
DATASET="https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.8v212"
ALLOWED_HOSTS={"datadryad.org"}
MAX_BYTES=2_500_000
UA="ALAN-CC0-moth-feeding-source-audit/1.0 (https://github.com/zuizui0223/alan)"
ENDPOINTS=[
    DATASET,
    "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.8v212/versions",
    "https://datadryad.org/downloads/file_stream/27006",
    "https://datadryad.org/api/v2/files/27006/download"
]
EXPECTED={"Day","Block","Treatment","Species","Sex","unique.id","Obs","Size"}


def fetch(url:str)->tuple[dict,bytes]:
    record={"http_status":None,"size_bytes":None,"sha256":None,
            "error":None,"requested_url":url.split("?")[0]}
    try:
        if urllib.parse.urlsplit(url).hostname not in ALLOWED_HOSTS:
            record["error"]="HOST_NOT_AUTHORIZED"
            return record,b""
        req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json,text/plain,text/csv,*/*"})
        with urllib.request.urlopen(req,timeout=32) as response:
            data=response.read(MAX_BYTES+1)
            record["http_status"]=response.status
            final=response.geturl()
            if urllib.parse.urlsplit(final).hostname not in ALLOWED_HOSTS:
                record["error"]="REDIRECT_HOST_NOT_DRYAD"
                return record,b""
            record["content_type"]=response.headers.get("Content-Type")
        if len(data)>MAX_BYTES:
            record["error"]="OVER_2P5MB_CAP"
            return record,b""
        record["size_bytes"]=len(data)
        record["sha256"]=hashlib.sha256(data).hexdigest()
        return record,data
    except urllib.error.HTTPError as err:
        record["http_status"]=err.code
        record["error"]="HTTP_"+str(err.code)
    except (urllib.error.URLError,ValueError,TimeoutError) as err:
        record["error"]=type(err).__name__+":"+str(err)[:120]
    return record,b""


def file_links(j:dict)->list[str]:
    links=[]
    for f in (j.get("_embedded",{}).get("stash:files",[]) if isinstance(j.get("_embedded"),dict) else []):
        for v in (f.get("_links",{}) or {}).values():
            if isinstance(v,dict) and v.get("href"):links.append(v["href"])
    for f in j.get("files",[]) or []:
        for key in ("downloadUrl","download_url","url"):
            if isinstance(f.get(key),str):links.append(f[key])
    return links


def parse_original(data:bytes)->tuple[list[dict],dict]:
    decoded=data.decode("utf-8-sig","replace")
    header=next((l for l in decoded.splitlines() if l.strip()),"")
    if not header:
        raise ValueError("Source data empty")
    likely=max(("\t",";",","),key=lambda t:header.count(t))
    if header.count(likely)==0:
        raise ValueError("Cannot safely recognize source column delimiter")
    reader=csv.DictReader(io.StringIO(decoded),delimiter=likely)
    names=set((reader.fieldnames or []))
    if not EXPECTED<=names:
        raise ValueError("Unexpected authentic Dryad schema; found "+str(sorted(names)))
    rows=list(reader)
    if not rows:
        raise ValueError("No authentic original data rows")
    return rows,{"rows_longitudinal":len(rows),"delimiter":repr(likely),"source_column_names":reader.fieldnames}


def derive(rows:list[dict])->dict:
    moth=collections.defaultdict(lambda:{"obs":0,"feeding":0,"species":set(),"treatment":set(),"blocks":set(),"days":set()})
    unknown=[]
    for row in rows:
        iid=str(row["unique.id"]).strip()
        if not iid:raise ValueError("Missing source moth individual ID")
        x=str(row["Obs"]).strip()
        if x in ("","NA","NaN","nan"):
            unknown.append(iid)
            continue
        if x not in ("0","1","0.0","1.0"):
            raise ValueError("Original Obs is not binary: "+repr(x))
        q=moth[iid]
        q["obs"]+=1;q["feeding"]+=int(float(x))
        q["species"].add(str(row["Species"]).strip())
        q["treatment"].add(str(row["Treatment"]).strip())
        q["blocks"].add(str(row["Block"]).strip())
        q["days"].add(str(row["Day"]).strip())
    byspecies=collections.defaultdict(list)
    contradictions=0
    for iid,q in moth.items():
        if len(q["species"])!=1 or len(q["treatment"])!=1 or len(q["blocks"])!=1:
            contradictions+=1
            continue
        taxon=next(iter(q["species"]))
        arm=next(iter(q["treatment"]))
        block=next(iter(q["blocks"]))
        if q["obs"]>0:
            byspecies[(taxon,arm,block)].append(q["feeding"]/q["obs"])
    if contradictions:
        raise ValueError("Source moth ID occurs in multiple experimental species/lights/blocks; need original identity audit: "+str(contradictions))
    alltaxa=sorted(set(k[0] for k in byspecies))
    allar=sorted(set(k[1] for k in byspecies))
    dark=[a for a in allar if a.casefold() in ("dark","control","no light","unlit")]
    if len(dark)!=1:raise ValueError("Cannot identify one original dark arm: "+repr(allar))
    dark=dark[0]
    means=collections.defaultdict(lambda:collections.defaultdict(dict))
    for (taxon,arm,block),values in byspecies.items():
        means[taxon][arm][block]=statistics.mean(values)
    results=[]
    for taxon in alltaxa:
        for arm in allar:
            samples=[v for (s,t,b),vs in byspecies.items() if s==taxon and t==arm for v in vs]
            paired=[]
            if arm!=dark:
                for block,x in means[taxon][arm].items():
                    control=means[taxon][dark].get(block)
                    if control is not None:paired.append(x-control)
            results.append({
                "taxon":taxon,"lamp_treatment":arm,
                "n_unique_moths":len(samples),
                "feeding_fraction_per_moth_mean":None if not samples else statistics.mean(samples),
                "n_blocks_observed":len(means[taxon][arm]),
                "blocks_with_same_species_dark_pair":len(paired),
                "matched_block_mean_difference_vs_dark":None if not paired else statistics.mean(paired),
                "leave_one_paired_block_min_max":None if len(paired)<2 else [
                    min(statistics.mean([q for k,q in enumerate(paired) if k!=i]) for i in range(len(paired))),
                    max(statistics.mean([q for k,q in enumerate(paired) if k!=i]) for i in range(len(paired)))],
                "matched_block_inferential_p":None,
                "plant_reproductive_measurement":False
            })
    return {
      "total_original_longitudinal_observations":len(rows),
      "source_missing_binary_obs":len(unknown),
      "unique_individual_moths_with_binary_data":len(moth),
      "moth_species_n":len(alltaxa),
      "moth_species_names":alltaxa,
      "lighting_treatments":allar,
      "original_dark_label":dark,
      "individuals_not_records_are_independent_observational_units":True,
      "results":results,
      "original_lamp_assignment_cage_randomization_rechecked":False,
      "confirmed_light_sensitivity_of_pollination_efficiency":False,
      "confirmed_reproductive_meta_beta_V":None
    }


def run()->dict:
    meta={
      "schema":"ALAN_MothFeeding2017_Dryad_public_source_gate_v1",
      "generated_utc":datetime.now(timezone.utc).isoformat(),
      "source_article_doi":"10.1098/rsbl.2016.0874",
      "dataset_doi":DOI,
      "expected_file":"Data_Moth_Foraging.txt",
      "licence":"Source Dryad page public CC0",
      "source_receipts":[],
      "actual_original_source_data_acquired":False,
      "original_data_analysis":None,
      "biological_scope":"Feeding frequency under ALAN, NOT effective pollen deposition or seed/fruit response",
      "meta_ALAN_reproduction_slope":None
    }
    possibilities=ENDPOINTS.copy()
    visited=set()
    while possibilities and len(visited)<14:
        url=possibilities.pop(0)
        if url in visited:continue
        visited.add(url)
        info,data=fetch(url)
        meta["source_receipts"].append(info)
        if not data:continue
        if data.lstrip()[:1]==b"{":
            try:
                j=json.loads(data)
                if isinstance(j,dict):
                    for u in file_links(j):
                        u=urllib.parse.urljoin(url,u)
                        if urllib.parse.urlsplit(u).hostname in ALLOWED_HOSTS and u not in visited:
                            possibilities.append(u)
                    if j.get("id") and "versions" in url:
                        for version in j.get("_embedded",{}).get("stash:versions",[]) or []:
                            if isinstance(version,dict) and version.get("id"):
                                possibilities.append("https://datadryad.org/api/v2/versions/"+str(version["id"])+"/files")
                    if "versionNumber" in j and j.get("id"):
                        possibilities.append("https://datadryad.org/api/v2/versions/"+str(j["id"])+"/files")
                    if "fileSize" in j and j.get("id"):
                        possibilities.append("https://datadryad.org/api/v2/files/"+str(j["id"])+"/download")
            except (ValueError,TypeError):
                pass
            continue
        # Safe title+schema gate: no HTML, ZIP or PDF mistaken for longform.
        if b"<html" in data[:500].lower() or data.startswith(b"%PDF") or data.startswith(b"PK"):
            continue
        try:
            rows,structure=parse_original(data)
            ans=derive(rows)
            meta.update(actual_original_source_data_acquired=True,
                        original_source_file_sha256=info["sha256"],
                        original_data_structure=structure,
                        original_data_analysis=ans)
            break
        except ValueError as err:
            info["parse_rejection"]=str(err)[:250]
    destination=pathlib.Path("results/alan_mothfeeding2017_dryad_individual_light_20261115.json")
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
      "original_source_acquired":meta["actual_original_source_data_acquired"],
      "checked_urls":len(visited),
      "statuses":[{"status":r["http_status"],"path":r["requested_url"][-100:]} for r in meta["source_receipts"]],
      "analysis":meta["original_data_analysis"]
    },ensure_ascii=False))
    return meta

if __name__=="__main__":
    run()
