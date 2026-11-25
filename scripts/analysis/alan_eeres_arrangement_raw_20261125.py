"""Source-gated actual arrangement analysis, pending original workbook.
Usage: python scripts/analysis/alan_eeres_arrangement_raw_20261125.py DISCONNECTED_alldata.xlsx --out result.json
The public Dryad DOI is 10.5061/dryad.rfj6q57k5. This file is not included here.
Do not use a fitted upstream mean as the missing control branch.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t

SHA="496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f"
def valid_taxa(df,offset,n):
    cols=list(df.columns[offset:])
    if len(cols)!=n or len(set(cols))!=n:
        raise ValueError("TAXON_SCHEMA_STOP")
    x=df[cols].apply(pd.to_numeric,errors="raise")
    if not np.isfinite(x.to_numpy()).all() or (x<0).any().any() or (x%1!=0).any().any():
        raise ValueError("TAXON_COUNT_STOP")
    return x

def require_sizes(df,keys,n,groups):
    s=df.groupby(keys,dropna=False).size()
    if len(s)!=groups or not s.eq(n).all():
        raise ValueError("PHYSICAL_GRAIN_STOP "+str(keys))

def calc(d,b):
    if len(d)!=144 or len(b)!=108:
        raise ValueError("SOURCE_ROWS_STOP")
    for df in (d,b):
        if df.ID.isna().any() or df.ID.duplicated().any():
            raise ValueError("SAMPLE_ID_STOP")
        for k in ("HT","flume"):
            df[k]=df[k].astype(str).str.strip().str.lower()
        if set(df.HT)!={f"ht{i}" for i in range(1,7)} or df.flume.nunique()!=18:
            raise ValueError("HEADER_FLUME_STOP")
    x=valid_taxa(d,8,38); y=valid_taxa(b,6,40)
    if int(x.to_numpy().sum())!=5564 or int(y.to_numpy().sum())!=11326:
        raise ValueError("PUBLISHED_TOTAL_STOP")
    if set(d.time)!={"T1","T2","T3","T4"}:
        raise ValueError("TIME_STOP")
    require_sizes(d,["HT","flume","time"],2,72)
    if d.duplicated(["HT","flume","time","side"]).any():
        raise ValueError("TRIBUTARY_DUPLICATE_STOP")
    expected={"CT":["CT","CT"],"SS":["F","L"],"CS":["CT","L+F"]}
    for (_,_,_),z in d.groupby(["HT","flume","time"]):
        arm=z.treatment_down.unique()
        if len(arm)!=1 or arm[0] not in expected or sorted(z.treatment_up)!=expected[arm[0]]:
            raise ValueError("PLACEMENT_STOP")
    mapping=d[["HT","flume","treatment_down"]].drop_duplicates()
    if len(mapping)!=18 or mapping.duplicated(["HT","flume"]).any():
        raise ValueError("TREATMENT_DRIFT_STOP")
    require_sizes(mapping,["HT","treatment_down"],1,18)
    if set(mapping.treatment_down)!=set(expected):raise ValueError("NETWORK_ARMS_STOP")
    # Reproduce authors' combined tributary count and union richness at each flume-week.
    dat=pd.concat([d[["HT","flume","time","treatment_down"]],x],axis=1)
    a=dat.groupby(["HT","flume","time","treatment_down"])[list(x.columns)].sum()
    events=pd.DataFrame({"abundance":a.sum(axis=1),"richness":(a>0).sum(axis=1)}).reset_index()
    require_sizes(events,["HT","flume","treatment_down"],4,18)
    drift=events.groupby(["HT","treatment_down"])[["abundance","richness"]].mean().reset_index().rename(columns={"treatment_down":"arm"})
    # Downstream benthos is a separate terminal sample, not drift.
    bd=b[b.section.astype(str).str.strip().str.lower().eq("down")].copy()
    if len(bd)!=36:raise ValueError("BENTHOS_DOWN_ROWS_STOP")
    bd["arm"]=bd.treatment.replace({"CT_down":"CT"})
    check=bd.merge(mapping,on=["HT","flume"],how="left",validate="many_to_one")
    if not check.arm.eq(check.treatment_down).all():raise ValueError("CROSS_SHEET_ARM_STOP")
    require_sizes(bd,["HT","flume","arm"],2,18)
    bd["abundance"]=y.loc[bd.index].sum(axis=1);bd["richness"]=(y.loc[bd.index]>0).sum(axis=1)
    benthos=bd.groupby(["HT","arm"])[["abundance","richness"]].mean().reset_index()
    outputs=[]
    for name,g in (("summed_tributary_drift",drift),("downstream_benthos",benthos)):
        for endpoint in ("abundance","richness"):
            p=g.pivot(index="HT",columns="arm",values=endpoint)
            if p.shape!=(6,3) or p.isna().any().any() or set(p.columns)!=set(expected):
                raise ValueError("BLOCK_ARM_SUPPORT_STOP")
            for num,den in (("CS","SS"),("SS","CT"),("CS","CT")):
                a=p[num].to_numpy();c=p[den].to_numpy()
                m=a.mean();n=c.mean()
                rec=dict(stage=name,endpoint=endpoint,numerator=num,denominator=den,n_blocks=6,mean_num=float(m),mean_den=float(n),block_num=a.tolist(),block_den=c.tolist())
                diff=a-c
                rec["paired_difference"]=float(diff.mean())
                if m>0 and n>0:
                    score=a/m-c/n; se=float(score.std(ddof=1)/np.sqrt(6)); theta=float(np.log(m/n))
                    rec.update(RR=float(np.exp(theta)),SE_logRR=se,CI95=np.exp(theta+np.array([-1,1])*t.ppf(.975,5)*se).tolist())
                else:rec["status"]="NONPOSITIVE_ARM_MEAN_NO_RR"
                outputs.append(rec)
    if abs(float(drift.abundance[drift.arm.eq("SS")].mean())-74.58333333333333)>1e-8:
        raise ValueError("FROZEN_SS_MEAN_STOP")
    return dict(status="SOURCE_REPLAYED_ARRANGEMENT_EXPLORATORY",source_sha256=SHA,independent_headers=6,results=outputs,limitation="Posthoc t5 paired-delta intervals, no family-wise correction; combined tributary drift is not independently sampled downstream flux; downstream benthos terminal; same programme, no new meta k.")

def main():
    p=argparse.ArgumentParser();p.add_argument("xlsx",type=Path);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if hashlib.sha256(a.xlsx.read_bytes()).hexdigest()!=SHA:raise ValueError("ORIGINAL_SHA_STOP")
    d=pd.read_excel(a.xlsx,sheet_name="drift");b=pd.read_excel(a.xlsx,sheet_name="benthos")
    result=calc(d,b)
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
if __name__=="__main__":main()
