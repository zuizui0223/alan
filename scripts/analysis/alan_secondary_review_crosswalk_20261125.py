"""Hash-gated review bibliography intake, not effect-size adoption.
Usage: python scripts/analysis/alan_secondary_review_crosswalk_20261125.py
       --source seymoure|czarnecka|czarnecka_s2 --input FILE --out DIR
       [--fetch instead of --input for the two previously verified file IDs]
Do not use a Data DOI or a review's study number as independent primary study k.
"""
import argparse, csv, hashlib, io, json, re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from openpyxl import load_workbook

SOURCES={
 "czarnecka":{"doi":"10.5061/dryad.xwdbrv1q0","file":"Data_S1_Czarnecka_et_al..xlsx","url":"https://datadryad.org/downloads/file_stream/4426061","sha":"8fae8a0244047ab057808ea321962bf848cff42b6fd521f34b37a360efa154d6"},
 "seymoure":{"doi":"10.6084/m9.figshare.c.6837630","file":"rstb20220356_si_002.xlsx","url":"https://figshare.com/ndownloader/files/42366918","sha":"a4c105c50f7b53cb2c58b287d065b79432b6ec7c308c35b1c4b873e78a4a359e"},
 "czarnecka_s2":{"doi":"10.1111/gcb.70569","file":"gcb70569-sup-0002-DataS2.xlsx","url":None,"sha":None}
}
ROUTING=Path("database/alan_experiments/unified_meta_intake_20261125/report_routing.csv")
HEADERS={"doi","doi link","article doi","paper doi","study doi","publication doi","reference doi","doi of article","doi of study"}
def norm(s):
    s=str(s or "").strip().lower()
    s=re.sub(r"^https?://(?:dx\.)?doi\.org/","",s)
    s=re.sub(r"^doi:\s*","",s).rstrip(" .,;")
    return s if re.fullmatch(r"10\.\d{4,9}/[^\s<>]+",s) else ""
def label(s):
    return re.sub(r"\s+"," ",re.sub(r"[_-]+"," ",str(s or "").strip().lower()))
def title(s):
    return re.sub(r"[^\w]+"," ",str(s or "").casefold()).strip()
def load_routes(path):
    a=list(csv.DictReader(path.open(encoding="utf8",newline="")))
    if len(a)!=732 or len({r["report_key"] for r in a})!=732:raise ValueError("ROUTING_SNAPSHOT_STOP")
    return a
def get(source,path,fetch):
    if path:
        return path.read_bytes()
    if not fetch or not source["url"] or not source["sha"]:raise ValueError("LOCAL_FILE_REQUIRED")
    # Bounded official public retrieval, no credentials or access circumvention.
    with urlopen(Request(source["url"],headers={"User-Agent":"ALAN-source-audit/1.0"}),timeout=35) as r:
        b=r.read(8_000_001)
    if len(b)>8_000_000:raise ValueError("SIZE_STOP")
    return b
def run(data,kind,routes):
    src=SOURCES[kind];digest=hashlib.sha256(data).hexdigest()
    if src["sha"] and digest!=src["sha"]:raise ValueError("SHA256_STOP "+digest)
    if not data.startswith(b"PK\x03\x04"):raise ValueError("NOT_XLSX_STOP")
    wb=load_workbook(io.BytesIO(data),read_only=True,data_only=True)
    lookup={norm(x["report_key"]):x for x in routes if x["report_key"].startswith("doi:") and norm(x["report_key"])}
    bytitle=defaultdict(list)
    for x in routes:
        if title(x["title"]):bytitle[title(x["title"])].append(x["report_key"])
    out=[];schema=[];seymoure_keys=set();nscientific=0
    for ws in wb.worksheets:
        if kind=="seymoure" and ws.title!="database":continue
        rows=list(ws.iter_rows(values_only=True))
        # Only use clearly labelled publication DOI or title fields, never scan arbitrary effect cells.
        candidates=[(i,r) for i,r in enumerate(rows[:15]) if any(label(x) in HEADERS for x in r)]
        tcandidates=[(i,r) for i,r in enumerate(rows[:15]) if any(label(x) in {"article title","paper title","publication title","study title","title"} for x in r)]
        if kind=="seymoure":
            candidates=[(i,r) for i,r in candidates if "interaction_type" in [str(x or "").strip() for x in r]]
        schema.append({"sheet":ws.title,"nonempty_rows":sum(any(v is not None and str(v).strip() for v in r) for r in rows),"first_nonempty_rows":[[str(v)[:150] if v is not None else "" for v in r] for r in rows[:3]],"candidate_doi_headers":[i+1 for i,r in candidates],"candidate_title_headers":[i+1 for i,r in tcandidates]})
        if len(candidates)==1:hi,hr=candidates[0]
        elif len(candidates)==0 and len(tcandidates)==1 and kind!="seymoure":hi,hr=tcandidates[0]
        else:continue
        headers=[str(x or "").strip() for x in hr];dc=[i for i,x in enumerate(headers) if label(x) in HEADERS];tc=[i for i,x in enumerate(headers) if label(x) in {"article title","paper title","publication title","study title","title"}]
        for i,r in enumerate(rows[hi+1:],hi+2):
            if not any(v is not None and str(v).strip() for v in r):continue
            vals={k:str(v).strip() if v is not None else "" for k,v in zip(headers,r) if k}
            if kind=="seymoure":
                if not vals.get("interaction_type"):
                    expected={"interaction strength (IS) reduction":"29","IS enhancement":"28","strength sum":"57","temporal mismatch induced":"7","spatial mismatch induced":"27"}
                    if any(vals.get(k)!=v for k,v in expected.items()) or vals.get("DOI") or vals.get("Article Title"):raise ValueError("FOOTER_STOP")
                    continue
                nscientific+=1
            ds=sorted({norm(r[j]) for j in dc if j<len(r) and norm(r[j])})
            tt=next((str(r[j]).strip() for j in tc if j<len(r) and r[j] is not None and str(r[j]).strip()),"")
            if kind=="seymoure":
                if len(ds)>1:raise ValueError("CONFLICTING_DOI_STOP")
                if ds:seymoure_keys.add("doi:"+ds[0])
                else:
                    a,b,c=vals.get("Authors",""),vals.get("Article Title",""),vals.get("Publication Year","")
                    if not all((a,b,c)):raise ValueError("SEYMOUR_ID_STOP")
                    seymoure_keys.add("fallback:"+repr((a.casefold(),b.casefold(),c)))
            matches=[lookup[d]["report_key"] for d in ds if d in lookup]
            state="exact_doi_match" if len(ds)==1 and len(matches)==1 else "new_doi_candidate" if len(ds)==1 else "conflicting_doi_hold" if len(ds)>1 else "no_publication_doi_hold"
            # An exact normalized title is a manual-review hint only. Never merge automatically.
            hints=bytitle.get(title(tt),[]) if tt else []
            out.append(dict(source=kind,sheet=ws.title,excel_row=i,publication_doi=";".join(ds),title=tt,interaction_type=vals.get("interaction_type",""),status=state,matched_report=";".join(matches),title_match_hints=";".join(hints),final_eligibility="not_adjudicated",source_hash=digest))
    if kind=="seymoure" and (nscientific!=81 or len(seymoure_keys)!=79):raise ValueError("79/81_SOURCE_COUNT_STOP")
    if not src["sha"]:
        # S2 lacks an independently recorded byte hash. Nothing is automatically admitted.
        status="LOCAL_S2_SCHEMA_AND_IDENTITY_HINTS_UNVERIFIED_BYTES"
    elif not out:status="VERIFIED_BYTES_SCHEMA_ONLY_HEADER_MAPPING_REQUIRED"
    else:status="VERIFIED_BYTES_DOI_ROUTING_ONLY"
    known={r["publication_doi"] for r in out if r["status"]=="exact_doi_match"}
    unknown={r["publication_doi"] for r in out if r["status"]=="new_doi_candidate"}
    return out,{"source":kind,"sha256":digest,"source_status":status,"sheet_schema":schema,"routing_rows":len(out),"distinct_matched_dois":len(known),"distinct_unmatched_doi_candidates":len(unknown),"row_status_counts":dict(Counter(r["status"] for r in out)),"seymoure_verified_report_count":len(seymoure_keys) if kind=="seymoure" else None,"new_meta_k":0,"final_eligibility":False}
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",choices=SOURCES,required=True);p.add_argument("--input",type=Path);p.add_argument("--fetch",action="store_true");p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    try:
        b=get(SOURCES[a.source],a.input,a.fetch);rows,receipt=run(b,a.source,load_routes(ROUTING))
        if rows:
            with (a.out/"report_crosswalk.csv").open("w",encoding="utf8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    except Exception as e:
        receipt={"source":a.source,"source_status":"SOURCE_STOP","error":type(e).__name__+": "+str(e),"new_meta_k":0,"final_eligibility":False}
    (a.out/"receipt.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in receipt.items() if k!="sheet_schema"},ensure_ascii=False))
if __name__=="__main__":main()
