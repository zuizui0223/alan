#!/usr/bin/env python3
"""Obtain actual source figure/supplement availability for P1×P3 direct studies.

Search ONLY original publisher / institutional / official repository endpoints.
Do not mistake HTTP success for accessible numerical factorial treatment rows.
Write self-contained source provenance and conservative published outcomes.
"""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,html.parser,json,re,urllib.error,urllib.request
from pathlib import Path
from urllib.parse import urljoin,urlsplit
SOURCES={
 "lu_JPE_publisher":"https://academic.oup.com/jpe/advance-article/doi/10.1093/jpe/rtag160/8731923",
 "lu_JPE_chinese_journal":"https://www.jpe.ac.cn/EN/10.1093/jpe/rtag160",
 "lu_JPE_DOI":"https://doi.org/10.1093/jpe/rtag160",
 "lu_JPE_accepted_pdf":"https://academic.oup.com/jpe/advance-article-pdf/doi/10.1093/jpe/rtag160/68754751/rtag160.pdf",
 "lu_JPE_crossref":"https://api.crossref.org/works/10.1093/jpe/rtag160",
 "lu_JPE_openalex":"https://api.openalex.org/works/https://doi.org/10.1093/jpe/rtag160",
 "mac_original_openaccess_pdf":"https://eprints.whiterose.ac.uk/id/eprint/149121/1/Macgregor_et_al_2019_Ecosphere.pdf",
 "mac_publisher_appendix_pdf":"https://esajournals.onlinelibrary.wiley.com/doi/supinfo/10.1002/ecs2.2550/ecs22550-sup-0001-AppendixS1.pdf",
 "mac_publisher_appendix_alt":"https://esajournals.onlinelibrary.wiley.com/doi/suppl/10.1002/ecs2.2550/supinfo/ecs22550-sup-0001-AppendixS1.pdf",
 "mac_original_code_repository":"https://raw.githubusercontent.com/CallumJMacgregor/BioassayPaper/master/FinalScripts.R",
}
CC_OPEN_PDFS={"mac_original_openaccess_pdf","lu_JPE_accepted_pdf"}
class Links(html.parser.HTMLParser):
 def __init__(self):
  super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  if tag.lower()!="a":return
  d=dict(attrs);u=d.get("href","")
  if u and any(k in (u+" "+d.get("title","")).lower() for k in
             ["rtag160","supplement","data","appendix","pdf","jpe"]):
   self.links.append({"href":u,"title":d.get("title","")})

def get(name,url):
 entry={"id":name,"requested_url":url}
 try:
  req=urllib.request.Request(url,headers={
   "User-Agent":"Mozilla/5.0 (compatible; ALAN-source-validation/1.0)",
   "Accept":"text/html,application/pdf,application/vnd.openxmlformats-officedocument.wordprocessingml.document,application/json,*/*"})
  with urllib.request.urlopen(req,timeout=30) as r:
   payload=r.read(6_500_001)
   entry.update(http_status=r.status,final_url=r.url,
     content_type=r.headers.get("content-type"),content_length_header=r.headers.get("content-length"))
  if len(payload)>6_500_000:raise ValueError("source exceeded 6.5MB bound")
  entry.update(bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),
     payload_format=("PDF" if payload.startswith(b"%PDF") else
     "DOCX_ZIP" if payload[:2]==b"PK" and "docx" in str(url).lower() else
     "JSON" if payload.lstrip().startswith((b"{",b"[")) else "HTML_OR_TEXT"))
  if entry["payload_format"]=="HTML_OR_TEXT":
   doc=payload.decode("utf-8",errors="replace")
   l=Links();l.feed(doc)
   entry["candidate_primary_links"]=[{**it,"absolute":urljoin(entry["final_url"],it["href"])}
      for it in l.links[:200]]
   entry["contains_seed"]=bool(re.search(r"seed\s*set|seedset|结籽",doc,re.I))
  if entry["payload_format"]=="JSON":
   j=json.loads(payload);entry["top_keys"]=list(j)[:25] if isinstance(j,dict) else []
   if name.endswith("crossref"):
    v=j.get("message",{})
    entry["metadata_links"]=v.get("link",[])[:20]
   if name.endswith("openalex"):
    entry["best_oa_location"]=j.get("best_oa_location")
    entry["locations"]=[{"pdf_url":v.get("pdf_url"),"landing_page_url":v.get("landing_page_url"),"source":v.get("source",{}).get("display_name")} for v in j.get("locations",[])[:15]]
  if entry["payload_format"]=="PDF" and name in CC_OPEN_PDFS:
   Path("source_probe").mkdir(exist_ok=True)
   (Path("source_probe")/(name+".pdf")).write_bytes(payload)
   entry["saved_author_paper_pdf"]=True
  entry["retrieval_status"]="VERIFIED_BYTES"
 except Exception as e:
  entry.update(retrieval_status="UNAVAILABLE",error_type=type(e).__name__,
               http_status=getattr(e,"code",None),error=str(e)[:350])
 return entry

def main(out):
 dest=Path(out)
 if dest.exists():raise FileExistsError("No overwrite")
 dest.parent.mkdir(parents=True,exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  results=list(pool.map(lambda it:get(*it),SOURCES.items()))
 response={"schema":"ALAN_Mac2019_Lu2026_P1P3_primary_fulltext_data_access_v1",
  "sources":results,"papers":[
  {"doi":"10.1002/ecs2.2550","species":"Silene latifolia",
   "experiment":{"runs":4,"source_flowers":417,"source_individual_female_plants":169,
     "four_pollinator_access_levels":True,
     "five_lighting_levels":True},
   "published_results":{"lighting_main_effect_pollination_success":"SIGNIFICANT",
     "caging_main_effect_pollination_success":"SIGNIFICANT",
     "LIGHT_X_CAGING_pollination_success":"NON_SIGNIFICANT_REPORTED_NO_P_SE",
     "lighting_main_effect_seed_count_full_model":"NON_SIGNIFICANT",
     "LIGHT_X_CAGING_seed_count_full_model":"NON_SIGNIFICANT_REPORTED"},
   "peer_reviewed":True,
   "effect_modifier_numerical_beta_or_SE_retrieved":False,
   "source":"published article DOI 10.1002/ecs2.2550, main text results pp.7,12 (PDF)" },
  {"doi":"10.1093/jpe/rtag160",
   "species":["Pedicularis kansuensis","Euphrasia regelii"],
   "published_abstract_alpine_manipulation":{
     "reported_Lepidoptera_visit_reduction_fraction":0.491,
     "reported_night_pollen_reduction_fraction":0.538,
     "reported_seed_set_reduction_fraction_Pedicularis":0.56,
     "reported_autonomous_selfing_seed_share_Euphrasia":0.65},
   "reported_Pedicularis_seed_set_lit_dark_ratio_from_percent":0.44,
   "reported_Pedicularis_seed_set_log_ratio_from_percent":__import__("math").log(.44),
   "same_species_light_x_diel_numeric_cell_means_retrieved":False,
   "Lu_experiment_sample_and_variance_retrieved":False,
   "source":"DOI 10.1093/jpe/rtag160 accepted manuscript abstract"},
  ],
  "direct_P1xP3_interaction_coefficients_newly_estimated":0,
  "authors_published_non_significant_interaction_is_NOT_equivalence":True,
  "original_number_of_citable_studies":2,
  "limitations":["No fictional access-by-lamp four/five-way cell means",
    "No model SE/p manufactured from article narrative",
    "Percent reduction is not a diel×ALAN contrast unless reference cell defined",
    "Published manuscript supplementary docx can contain statistical tables but not assumed raw observations",
    "A failed URL is not proof that author-held data do not exist"]}
 dest.write_text(json.dumps(response,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
 for r in results:print("P1_P3_PRIMARY_SOURCE_GATE",r["id"],r["retrieval_status"],r.get("http_status"),r.get("payload_format"),r.get("bytes"))
 print("P1_P3_DIRECT_INTERACTION_NOT_YET_REESTIMATED")
if __name__=="__main__":
 ap=argparse.ArgumentParser();ap.add_argument("--out",required=True)
 main(ap.parse_args().out)
