#!/usr/bin/env python3
"""Original-paper Grognuz 2026 16 flowering plant label × UK moisture support.

Publication: Grognuz, Gisler & Knop, Biological Conservation 314:111650,
DOI 10.1016/j.biocon.2025.111650 (publication date 2026).
Authors' figure 4 has sixteen labels; two genus-only labels, one accepted
species binomial with no recorded UK F, and one unresolved synonym. No
figure-pixel effect sizes are reinterpreted as source plot observations.

Original LIGHT study:
 2022: six lit/dark site pairs, but 1 dark site excluded from phenology
 2023: eight lit/dark site pairs
 Original PHENOLOGY plots: 108 (44 in 2022 + 64 in 2023).
 Experimental unit is lit/dark physical SITE, not plot×weekly×species.

Original paper's Mendeley dataset DOI is reserved for version TWO
(10.17632/xjmn5jhxyn.2), while publicly indexed web page currently
shows version ONE 10.17632/xjmn5jhxyn.1.
Audit both versions through official current api.data.mendeley.com public
API; a 401 / 403 / version conflict is a source STOP not an inferred zero.

UK traits come from RBG Kew's BIFloraExplorer (Git blob SHA1 pinned).
No author figure species alias is invented, no 'sp.' genus moisture.
The phenology effect can only be fit from true original site/plot data,
never by treating 16 species x 4 endpoints as independent ALAN studies.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,urllib.request,urllib.error,concurrent.futures
from pathlib import Path
SCHEMA="alan_Grognuz2026_Fig4_original16_KewF_Mendeley_access_support_v1"
PAPER_DOI="10.1016/j.biocon.2025.111650"
PUBPDF="https://ira.agroscope.ch/en-US/Page/Einzelpublikation/Download?einzelpublikationId=68011"
SOURCE_PUBLIC="10.17632/xjmn5jhxyn.1"
SOURCE_RESERVED="10.17632/xjmn5jhxyn.2"
K_GIT_BLOB="e6a1da929af81f370d2a1b79177b489e3e573a06"
K_RAW_URL="https://raw.githubusercontent.com/RBGKew/BIFloraExplorer/master/data-raw/BI_MAIN.csv"
M_ID="xjmn5jhxyn"
M_API="https://api.data.mendeley.com"
PLANTS=(
    "Achillea millefolium",
    "Centaurea cyanus",
    "Sinapis arvensis",
    "Daucus carota",
    "Papaver rhoeas",
    "Sonchus sp.",
    "Centaurea jacea",
    "Leucanthemum vulgare",
    "Malva moschata",
    "Silene pratensis",
    "Trifolium incarnatum",
    "Echium vulgare",
    "Epilobium sp.",
    "Veronica persica",
    "Matricaria chamomilla",
    "Trifolium pratense",
)
assert len(PLANTS)==16 and len(set(PLANTS))==16
def check(ok,msg):
    if not ok:raise RuntimeError("GROGNUZ_ORIGINAL16_SOURCE_GATE_STOP: "+msg)
def git_blob(raw):
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()
def get(url,maxsize=4<<20,headers=None):
    h={"User-Agent":"ALAN-research-Grognuz-original-source-audit/2026",
       "Accept":"application/json,application/vnd.mendeley-public-dataset.1+json,*/*"}
    if headers:h.update(headers)
    req=urllib.request.Request(url,headers=h)
    with urllib.request.urlopen(req,timeout=28) as response:
        blob=response.read(maxsize+1)
        if len(blob)>maxsize:raise ValueError("SOURCE_RESPONSE_SIZE_LIMIT_STOP")
        return blob,response.geturl(),response.headers.get("Content-Type","")
def kew_traits(blob):
    check(git_blob(blob)==K_GIT_BLOB,"original RBG Kew Git blob bytes changed")
    reader=csv.DictReader(io.StringIO(blob.decode("latin1")))
    required={"taxon_name_binom","F_PLANTATT","F_Doring",
              "taxonomic_status","aggregate"}
    check(required.issubset(reader.fieldnames or []),"original Kew F schema changed")
    target=set(PLANTS)
    results={};rows=0
    for line in reader:
        rows+=1
        species=line["taxon_name_binom"].strip()
        if species not in target:continue
        def parse(x):
            x=x.strip()
            if x in ("","NA","NaN","x","X"):return None
            num=float(x)
            check(1<=num<=12,"original F score out of range")
            return num
        f=parse(line["F_PLANTATT"])
        d=parse(line["F_Doring"])
        if species in results:
            check(results[species]["F_PLANTATT"]==f and results[species]["F_Doring"]==d,
                "duplicate conflicting Kew trait")
        results[species]={
          "F_PLANTATT":f,"F_Doring":d,
          "original_Kew_taxonomic_status":line["taxonomic_status"].strip(),
          "original_Kew_aggregate":line["aggregate"].strip()}
    check(rows>=3000,"Kew full flora unexpectedly incomplete")
    records=[]
    for name in PLANTS:
        val=results.get(name)
        ambiguous=name.endswith("sp.")
        records.append({"paper_Figure4_label":name,
             "Kew_exact_original_binomial":name if val else None,
             "not_a_species_level_label":ambiguous,
             "Kew_F_PLANTATT":val["F_PLANTATT"] if val else None,
             "Kew_F_Doring":val["F_Doring"] if val else None,
             "Kew_taxonomic_status":val["original_Kew_taxonomic_status"] if val else None,
             "source_Kew_missing_or_taxon_unresolved":val is None or val["F_PLANTATT"] is None,
             "no_genus_or_synonym_trait_assigned":True})
    return {"source_original_Kew_rows":rows,
            "original_Kew_source_Git_blob_SHA1":K_GIT_BLOB,
            "source_published_fig4_labels":len(records),
            "figure4_source_name_and_exact_moisture":records,
            "Kew_exact_binomial_found":sum(x["Kew_exact_original_binomial"] is not None for x in records),
            "F_PLANTATT_scored":sum(x["Kew_F_PLANTATT"] is not None for x in records),
            "F_Doring_scored":sum(x["Kew_F_Doring"] is not None for x in records),
            "F_PLANTATT_scored_wet_gt7":sum(x["Kew_F_PLANTATT"] is not None and
                                               x["Kew_F_PLANTATT"]>7 for x in records),
            "F_Doring_scored_wet_gt7":sum(x["Kew_F_Doring"] is not None and
                                             x["Kew_F_Doring"]>7 for x in records),
            "F_PLANTATT_distinct_scores":sorted({x["Kew_F_PLANTATT"] for x in records if x["Kew_F_PLANTATT"] is not None}),
            "F_Doring_distinct_scores":sorted({x["Kew_F_Doring"] for x in records if x["Kew_F_Doring"] is not None}),
            "Kew_F_are_species_habitat_indicators_not_actual_site_hydrology":True}
def probe(label,url):
    try:
        blob,final,ctype=get(url)
        result={"label":label,"url":url,"status":"HTTP_OK",
                "bytes":len(blob),"final_source_url":final,
                "http_content_type":ctype,
                "response_SHA256":hashlib.sha256(blob).hexdigest()}
        if "json" in ctype.lower() or blob.lstrip().startswith((b"{",b"[")):
            content=json.loads(blob)
            result["format"]="source_json"
            if isinstance(content,dict):
                result["top_level_fields"]=list(content.keys())[:50]
                result["source_name"]=content.get("name")
                result["version"]=content.get("version")
                result["id"]=content.get("id")
                doi=content.get("doi",{})
                result["doi"]=doi.get("id") if isinstance(doi,dict) else doi
                files=content.get("files") or content.get("results") or content.get("data")
                if isinstance(files,list):
                    result["original_file_records"]=len(files)
                    result["file_metadata"]=[{
                        "filename":f.get("filename"),
                        "id":f.get("id"),
                        "size":(f.get("content_details") or {}).get("size",f.get("size")),
                        "sha256":(f.get("content_details") or {}).get("sha256_hash"),
                        "download_url_explicitly_supplied":bool(
                            (f.get("content_details") or {}).get("download_url"))
                    } for f in files[:80] if isinstance(f,dict)]
            else:
                result["source_array_records"]=len(content) if isinstance(content,list) else None
                if isinstance(content,list):
                    result["file_metadata"]=[{
                        "filename":f.get("filename"),
                        "id":f.get("id"),
                        "size":(f.get("content_details") or {}).get("size",f.get("size")),
                        "sha256":(f.get("content_details") or {}).get("sha256_hash"),
                        "download_url_explicitly_supplied":bool(
                            (f.get("content_details") or {}).get("download_url"))
                    } for f in content[:80] if isinstance(f,dict)]
        else:
            result["format"]="not_json_original_data_unverified"
        return result
    except Exception as e:
        return {"label":label,"url":url,"status":"SOURCE_ACCESS_STOP",
                "error_type":type(e).__name__,
                "http_code":getattr(e,"code",None)}
def mendeley_audit():
    endpoints=[
      ("dataset_v1",f"{M_API}/datasets/{M_ID}?version=1"),
      ("dataset_v2",f"{M_API}/datasets/{M_ID}?version=2"),
      ("latest_dataset",f"{M_API}/datasets/{M_ID}"),
      ("files_v1",f"{M_API}/datasets/publics/{M_ID}/files?version=1"),
      ("files_v2",f"{M_API}/datasets/publics/{M_ID}/files?version=2"),
      ("dataset_versions",f"{M_API}/datasets/{M_ID}/versions"),
      ("repository_v1_public_page",f"https://data.mendeley.com/datasets/{M_ID}/1"),
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures=[executor.submit(probe,label,url) for label,url in endpoints]
        results=[f.result() for f in futures]
    candidate=[x for x in results if x["status"]=="HTTP_OK" and
               x.get("format")=="source_json"]
    original_file_metadata=[r for r in candidate if r.get("file_metadata")]
    return {
        "source_article_reserved_DOI_v2":SOURCE_RESERVED,
        "source_public_webpage_DOI_v1":SOURCE_PUBLIC,
        "mendeley_authorized_API_does_not_imply_download_availability":True,
        "original_endpoints":results,
        "source_json_metadata_readable":len(candidate)>0,
        "verified_public_file_manifest_available":len(original_file_metadata)>0,
        "source_raw_experimental_phenology_rows_downloaded":False,
        "plant_by_site_original_light_response_fitted":False,
        "status":"PUBLIC_DATA_FILE_METADATA_SOURCE_GATE_PASS_NEEDS_RAW_FILES" if original_file_metadata else
                 "ORIGINAL_MENDELEY_PUBLISHED_FILE_SOURCE_STOP",
    }
def run():
    kew_raw,_,_=get(K_RAW_URL,maxsize=2<<20)
    K=kew_traits(kew_raw)
    M=mendeley_audit()
    check(K["source_published_fig4_labels"]==16,"Figure4 list changed")
    check(K["F_PLANTATT_scored"]==12,"original Kew exact indicator coverage changed")
    check(K["F_PLANTATT_scored_wet_gt7"]==0,"original source wet taxa unexpectedly present")
    return {"schema":SCHEMA,"publication_doi":PAPER_DOI,
            "article_source_Agroscope_PDF":PUBPDF,
            "original_phenology_2022_paired_light_sites":6,
            "original_phenology_2023_paired_light_sites":8,
            "original_phenology_2022_observed_sites":11,
            "original_phenology_2023_observed_sites":16,
            "source_phenology_plots_total":108,
            "original_2022_one_dark_site_missing_from_flowering":True,
            "species_and_F_source":K,"Mendeley_original_file_source_audit":M,
            "wet_vs_drier_species_comparison_supported":False,
            "original_ALAN_plant_phenology_moisture_moderator_computed":False,
            "source_reported_16_species_12_with_shift_original_author_result_only":True,
            "not_independent_16_species_x_4_phenology_endpoint_light_experiments":True,
            "no_plant_seed_set_measured":True}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    result=run()
    f=Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("GROGNUZ2026_FIG4_KewF_and_MENDELEY_SOURCE_GATE "+
       json.dumps({"F_coverage":{k:v for k,v in result["species_and_F_source"].items()
                                 if k!="figure4_source_name_and_exact_moisture"},
                   "original_published_plant_label_support":
                        result["species_and_F_source"]["figure4_source_name_and_exact_moisture"],
                   "Mendeley":result["Mendeley_original_file_source_audit"]},
                   sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
