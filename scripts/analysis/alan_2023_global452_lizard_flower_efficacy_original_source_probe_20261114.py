#!/usr/bin/env python3
"""Probe original GLOBAL 452 lizard–flower interactions and separate 2,838-lizard phylogeny.

Review: Justicia Correcher et al. 2023, Ecography DOI 10.1111/ecog.06425,
162 studies / 452 interaction pairs, 91 islands, 41 archipelagos, 31 mainland
countries, 173 lizard/311 plant species, 90 verified pollination interactions.
The original XLS Supporting Information (129.5 kB) is essential for *NEW*
study-clustered island/continent functional-role estimates. Without the XLS,
reported source numbers remain purely published results, no new regression.

Separate phylogenetic related dataset: 2023 Molecular Phylogenetics and
Evolution DOI 10.1016/j.ympev.2023.107839 and figshare 23460260, 2,838
lizard species with flower-visit and seed-dispersal positive markers.
DO NOT infer individual lizard species with absent record are biologically
incapable; record a historical source-positive/documentation coding.

No raw individual rare species locality rows or article text are republished.
"""
from __future__ import annotations
import argparse,hashlib,io,json,urllib.error,urllib.request
from pathlib import Path
SOURCES={
 "wiley_official_2023_S1": [
  "https://nsojournals.onlinelibrary.wiley.com/doi/suppl/10.1111/ecog.06425/supinfo/ecog12967-sup-0001-AppendixS1.xls",
  "https://onlinelibrary.wiley.com/doi/suppl/10.1111/ecog.06425/supinfo/ecog12967-sup-0001-AppendixS1.xls",
  "https://onlinelibrary.wiley.com/doi/supinfo/10.1111/ecog.06425/supinfo/ecog12967-sup-0001-AppendixS1.xls",
  "https://nsojournals.onlinelibrary.wiley.com/doi/supinfo/10.1111/ecog.06425/supinfo/ecog12967-sup-0001-AppendixS1.xls"
 ],
 "figshare_lizard_macro_2838": [
   "https://api.figshare.com/v2/articles/23460260",
   "https://api.figshare.com/v2/file/download/23460260"
 ]
}
def request(u,max_bytes=5000000,timeout=18):
    try:
        req=urllib.request.Request(u,headers={
             "User-Agent":"ALAN-public-original-island-lizard-interactions-research/1.0",
             "Accept":"application/vnd.ms-excel,application/json,*/*"})
        with urllib.request.urlopen(req,timeout=timeout) as f:
            b=f.read(max_bytes+1)
            if len(b)>max_bytes:
                return None,{"status":"source_byte_limit","max_bytes":max_bytes}
            return b,{"status":"retrieved_bytes","download_bytes":len(b),
                      "sha256":hashlib.sha256(b).hexdigest(),
                      "content_type":str(f.headers.get("Content-Type",""))[:90],
                      "final_hostname":urllib.request.urlparse(f.url).hostname
                         if hasattr(urllib.request,"urlparse") else None}
    except urllib.error.HTTPError as e:
        return None,{"status":"http_"+str(e.code)}
    except (urllib.error.URLError,TimeoutError,OSError,ValueError) as e:
        return None,{"status":"source_unavailable","error_type":type(e).__name__}

def xls_sheets(raw):
    import xlrd
    book=xlrd.open_workbook(file_contents=raw,logfile=io.StringIO())
    info=[]
    for sh in book.sheets():
        info.append({"name":sh.name,"nrows":sh.nrows,"ncols":sh.ncols,
            "first_row_columns":[str(x)[:125] for x in sh.row_values(0)[:70]] if sh.nrows else [],
            "data_rows_not_exported":True})
    return {"source_excel_parse":"ok","sheets":info}

def main(out):
    out=Path(out)
    if out.exists():raise ValueError("Existing original research receipt cannot be overwritten")
    w=[]
    recovered=None
    for url in SOURCES["wiley_official_2023_S1"]:
        raw,r=request(url)
        rec={"source_url":url,**r}
        if raw is not None and raw[:8]==b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1":
            try:
                m=xls_sheets(raw)
                rec.update(m)
                recovered={"source_url":url,"digest_sha256":r["sha256"],
                           "original_XLS_size":len(raw),"sheets":m["sheets"]}
            except Exception as e:
                rec["parse_error_type"]=type(e).__name__
                rec["parse_error_string"]=str(e)[:180]
            w.append(rec)
            if recovered:break
        else:
            if raw:
                rec["file_signature_first_8_hex"]=raw[:8].hex()
                rec["original_binary_XLS_confirmed"]=False
            w.append(rec)
    fig={}
    raw,rec=request(SOURCES["figshare_lizard_macro_2838"][0])
    fig["api_article"]=SOURCES["figshare_lizard_macro_2838"][0]
    fig["fetch_receipt"]=rec
    if raw:
        try:
            a=json.loads(raw)
            fig["title"]=str(a.get("title",""))[:240]
            fig["published_date"]=a.get("published_date")
            fig["version"]=a.get("version")
            fig["public_license"]=a.get("license")
            fig["files"]=[{k:v.get(k) for k in ("id","name","size","download_url",
                         "is_link_only","computed_md5") if k in v}
                         for v in a.get("files",[])]
            fig["original_file_list_count"]=len(fig["files"])
        except ValueError:
            fig["source_json_parse_error"]=True
    status=("ORIGINAL_GLOBAL_LIZARD_REVIEW_XLS_VERIFIED" if recovered else
            "ORIGINAL_REVIEW_XLS_UNAVAILABLE_OR_NOT_XLS_ONLY_PUBLISHED_GLOBAL_SUMMARY")
    result={
     "schema":"ALAN_lizard_global_2023_452_interactions_original_source_20261114_v1",
     "status":status,"source_review_DOI":"10.1111/ecog.06425",
     "reported_review_source_counts_NOT_REANALYZED":{
         "unique_primary_papers":162,"interaction_records":452,
         "lizard_species":173,"plant_species":311,
         "independent_physical_islands":91,"archipelagos":41,
         "mainland_countries":31,
         "records_reported_on_islands_percent_approximately":79,
         "records_with_effective_pollination_evidence":90,
         "effective_pollination_evidence_approx_percent":20,
         "species_with_documented_obligate_lizard_pollination":3},
     "official_original_supporting_filename":"ecog12967-sup-0001-AppendixS1.xls",
     "official_original_publisher_requests":w,
     "original_review_XLS_verified_metadata":recovered,
     "separate_phylogenetic_figshare_source":fig,
     "source_interpretation_gates":{
       "none_of_452_record_interactions_is_randomly_sampled_from_lizard_abundance":True,
       "plant_lizard_relationship_not_a_measured_ALAN_treatment":True,
       "positive_same_species_plant_pollination_record_not_local_alternative_at_night":True,
       "pollination_confirmed_in_original_review_not_all_independent_studies":True,
       "original_data_unauthorized_not_silently_substituted":True,
       "island_79pct_documentation_may_be_recording_bias_not_ecological_prevalence":True,
       "no_effective_light_by_island_lizard_interaction_estimated":True,
       "original_raw_plant_locality_record_rows_not_exported":True
     }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_LIZARD_ISLAND_REAL_GLOBAL_SOURCE_FILE_PROBE",json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();main(a.out)
