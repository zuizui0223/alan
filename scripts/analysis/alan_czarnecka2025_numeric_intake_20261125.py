"""Read authenticated Czarnecka Dryad source; reproduce numerical gates and study metadata crosswalk.
Usage: python scripts/analysis/alan_czarnecka2025_numeric_intake_20261125.py <original ZIP> --out <directory>
Dependency: openpyxl. No workbook modifications or source raw-file redistribution.
All output rows are candidates, not automatically admissible.
"""
import argparse, csv, hashlib, io, json, math, re, zipfile
from collections import defaultdict
from pathlib import Path
from openpyxl import load_workbook
SHA="8fae8a0244047ab057808ea321962bf848cff42b6fd521f34b37a360efa154d6"
README="f4f4e7005be4bcdb84baae13dc652d9ecd27b8cd5f8556a1251207d9b7d1411b"
DIRECT={"Consumption","Predation success","Parasitism","Pollination effectivity","Decomposition","Interspecific interactions","Predation risk"}
T=2.7764451051977987
def mean(a):return sum(a)/len(a)
def finite(x):return isinstance(x,(float,int)) and math.isfinite(x)
def doi(x):return str(x or "").strip().lower().rstrip(" .,;").replace(" ","")
def writecsv(path,rows):
    if not rows:return
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
def read_source(path):
    with zipfile.ZipFile(path) as z:
        b=z.read("Data_S1_Czarnecka_et_al..xlsx")
        assert hashlib.sha256(b).hexdigest()==SHA
        assert hashlib.sha256(z.read("README.md")).hexdigest()==README
    w=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    assert w.sheetnames==["Arkusz1"]
    it=w["Arkusz1"].values;hh=next(it)
    cols=[(i,str(v).strip()) for i,v in enumerate(hh) if v is not None]
    a=[dict(source_excel_row=i,**{k:r[j] for j,k in cols}) for i,r in enumerate(it,2) if any(v is not None for v in r)]
    w.close()
    assert len(a)==4905 and len({r["Study"] for r in a})==314 and len({r["ID"] for r in a})==4900
    return a
def analyze(rows):
    field=[r for r in rows if r["Study.Type"]=="Field"]
    assert len(field)==1676 and len({r["Study"] for r in field})==128
    for r in field:
        r["direct_category_route"]=r["Sub.Category1"] in DIRECT or r["Sub.Category"]=="Pollination"
        r["vascular_route"]=r["Phylum"]=="Tracheophyta"
        r["explicit_route"]=r["Study"] in {"068","138"} or (r["Study"]=="021" and str(r["Measured.effect"]).strip().lower()=="no of seeds per flower")
        r["focus"]=any(r[k] for k in ("direct_category_route","vascular_route","explicit_route"))
        m,s,n,M,S,N=(r[k] for k in ("Control_Mean","Control_SD","Control_N","Experimental_Mean","Experimental_SD","Experimental_N"))
        g=all(map(finite,(m,s,n,M,S,N))) and m>0 and M>0 and s>=0 and S>=0 and n>=2 and N>=2
        r["positive_mean_sd_n_gate"]=g
        r["nonnegative_support_sd_inconsistent_if_raw"]=g and (s>math.sqrt(n)*m+1e-9 or S>math.sqrt(N)*M+1e-9)
        r["summary_lnRR_candidate"]=math.log(M/m) if g else None
        r["final_eligibility"]="not_adjudicated"
    focus=[r for r in field if r["focus"]];good=[r for r in field if r["positive_mean_sd_n_gate"]];bad=[r for r in field if r["nonnegative_support_sd_inconsistent_if_raw"]]
    assert len(focus)==365 and len({r["Study"] for r in focus})==34 and len(good)==1486 and len(bad)==79
    bounds=[]
    for r in field:
        k=r["Study"];m=str(r["Measured.effect"])
        U=1 if (k=="021" and m.lower().strip().startswith("proportion of")) or k=="231" else 100 if k=="065" and m=="pollination success" else None
        if U is None:continue
        q={"ID":r["ID"],"study":k,"assumed_range":[0,U]}
        for a,c,s,n in [("dark","Control_Mean","Control_SD","Control_N"),("lit","Experimental_Mean","Experimental_SD","Experimental_N")]:
            mu,sd,N=r[c],r[s],r[n]
            ub=math.sqrt(N/(N-1)*mu*(U-mu)) if 0<=mu<=U and N>1 else None
            q[a]={"mean":mu,"SD":sd,"N":N,"max_SD_if_raw":ub,"fails":ub is not None and sd>ub+1e-8}
        bounds.append(q)
    assert len(bounds)==14 and all(q["dark"]["fails"] and q["lit"]["fails"] for q in bounds)
    x=[r for r in field if r["Study"]=="188"]
    assert len(x)==8 and [r["ID"] for r in x]==[f"188.{i}" for i in range(1,9)]
    a=[]
    for r in x:
        m,M,ds,ls=(r[k] for k in ("Control_Mean","Experimental_Mean","Control_SE","Experimental_SE"))
        assert r["Control_N"]==r["Experimental_N"]==5
        assert math.isclose(r["Control_SD"],ds*math.sqrt(5)) and math.isclose(r["Experimental_SD"],ls*math.sqrt(5))
        a.append(dict(ID=r["ID"],y=math.log(M/m),v=(ds/m)**2+(ls/M)**2))
    g=[]
    for i,k in enumerate(("common_alien","rare_alien","common_native","rare_native")):
        y=a[i+4]["y"]-a[i]["y"];v=a[i+4]["v"]+a[i]["v"]
        g.append(dict(group=k,y=y,v=v,RR=math.exp(y),CI=[math.exp(y-T*math.sqrt(v)),math.exp(y+T*math.sqrt(v))]))
    y=mean([r["y"] for r in g]);v=mean([math.sqrt(r["v"]) for r in g])**2
    out={"field_rows":len(field),"focus_rows":len(focus),"positive_gate":len(good),"negative_support_inconsistent_if_raw":len(bad),"bounded_proportions":bounds,"liu_candidate":{"arms":a,"contrasts":g,"mean_RR":math.exp(y),"upper_variance":v,"approx_t4_CI":[math.exp(y-T*math.sqrt(v)),math.exp(y+T*math.sqrt(v))]},"new_admitted_meta_k":0}
    return field,focus,out
def crosswalk(rows,route):
    grouped=defaultdict(list);owners=defaultdict(set)
    for r in rows:
        grouped[r["Study"]].append(r)
        if doi(r["DOI"]):owners[doi(r["DOI"])].add(r["Study"])
    lookup={str(r["report_key"]).lower():r for r in route}
    a=[]
    for key,g in sorted(grouped.items()):
        ds=sorted({doi(r["DOI"]) for r in g if doi(r["DOI"])})
        hits=[d for d in ds if "doi:"+d in lookup]
        if not ds:s="missing_DOI_hold"
        elif len(ds)!=1:s="multiple_DOIs_hold"
        elif len(owners[ds[0]])>1:s="shared_DOI_across_study_IDs_hold"
        elif key in {"085","101","264"}:s="nonarticle_or_malformed_DOI_hold"
        else:s="exact_DOI_candidate_match" if hits else "unmatched_DOI_candidate"
        a.append(dict(study_id=key,source_rows=len(g),field_rows=sum(r["Study.Type"]=="Field" for r in g),source_dois="|".join(ds),source_titles="|".join(sorted({str(r["Title"]) for r in g})),candidate_matches="|".join(hits),status=s,final_eligibility="not_adjudicated"))
    return a
def main():
    p=argparse.ArgumentParser();p.add_argument("zip",type=Path);p.add_argument("--out",type=Path,required=True);p.add_argument("--route",type=Path,default=Path("database/alan_experiments/unified_meta_intake_20261125/report_routing.csv"));a=p.parse_args()
    rows=read_source(a.zip);field,focus,out=analyze(rows)
    a.out.mkdir(parents=True,exist_ok=True)
    writecsv(a.out/"all_field_candidates.csv",field);writecsv(a.out/"focus_candidates.csv",focus)
    route=list(csv.DictReader(a.route.open(encoding="utf8",newline="")));assert len(route)==732
    writecsv(a.out/"all314_study_crosswalk.csv",crosswalk(rows,route))
    (a.out/"numeric_audit.json").write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
if __name__=="__main__":main()
