"""Screen wet-cohort plant traits using matched geographic and survey strata."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import calendar
import json
import re
import time
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from statsmodels.stats.multitest import multipletests
from alan_rf_decline_traits import compare

repo=Path(__file__).resolve().parents[2]
base=repo/"results/all_species_forest"
out=base/"wet_drivers"
out.mkdir(exist_ok=True)
cache=repo.parent/"gbif_annual/pladias_driver_accounts"
cache.mkdir(exist_ok=True)
w=pd.read_csv(base/"wet_reproduction/species.csv").merge(pd.read_csv(base/"decline_traits/moisture_evidence.csv.gz")[["sp","site_region","island_band","site_band"]],on="sp",validate="1:1")
assert w.candidate.sum()==58 and len(w)==363
def account(r):
    path=cache/(r.sp+".json")
    if not path.exists():
        for attempt in range(4):
            try:
                res=requests.get(r.url,timeout=40);res.raise_for_status()
                break
            except requests.RequestException:
                if attempt==3:raise
                time.sleep(2**attempt)
        page=BeautifulSoup(res.content,"html.parser")
        assert r.account_name in [h.get_text(" ",strip=True) for h in page.select("h1")]
        fields={}
        for tag in page.select("span.font-italic"):
            key=tag.get_text(" ",strip=True)
            if key and key!=":" and tag.parent.b:
                fields[key]=" ".join(tag.parent.b.stripped_strings)
        path.write_text(json.dumps(fields,ensure_ascii=False),encoding="utf-8")
    return dict(json.loads(path.read_text(encoding="utf-8")),sp=r.sp)
with ThreadPoolExecutor(6) as pool:
    raw=pd.DataFrame(list(pool.map(account,w.dropna(subset=["url","account_name"]).itertuples()))).set_index("sp")
w=w.set_index("sp")
values={}
def field(k):return raw[k].reindex(w.index).str.replace(r"\s+"," ",regex=True).str.strip()
def numeric(k):
    return pd.to_numeric(field(k).str.extract(r"^([+−-]?\d+(?:\.\d+)?)\s*(?:%|$)",expand=False).str.replace("−","-",regex=False),errors="coerce")
for k in ["Whole-community disturbance frequency indicator value","Whole-community disturbance severity indicator value","Herb layer disturbance frequency indicator value","Herb layer disturbance severity indicator value","Ecological specialization index for all vegetation types","Index of colonization success (ICS)","Index of colonization potential (ICP)","Optimum successional age [years]","Number of narrow habitats in which the taxon has its optimum","Occurrence frequency in vegetation plots","Mean percentage cover in vegetation plots","Clonal index","Lateral spreading distance by clonal growth","Persistence of the clonal growth organ","Number of clonal offspring","Size of the belowground bud bank (root buds included)","Depth of the belowground bud bank (root buds included)","Life strategy (Pierce method, C-score)","Life strategy (Pierce method, S-score)","Life strategy (Pierce method, R-score)"]:
    values[k]=numeric(k)
for axis in ["Light","Temperature","Moisture","Reaction","Nutrient","Salinity"]:
    k=axis+" indicator value"
    values[k]=pd.to_numeric(field(k).str.extract(r"^(\d+(?:\.\d+)?)\s*(?:[–-]|$)",expand=False),errors="coerce")
months={m:i for i,m in enumerate(calendar.month_name) if m}
def flowering(v):
    if pd.isna(v):return (np.nan,np.nan)
    # Only simple contiguous named month ranges; omit ambiguous alternatives or winter wrap.
    match=re.fullmatch(r"([A-Za-z]+)(?:\s*[–-]\s*([A-Za-z]+))?",v)
    if not match:return (np.nan,np.nan)
    a=months.get(match[1]);b=months.get(match[2] or match[1])
    return (a,b-a+1) if a and b and b>=a else (np.nan,np.nan)
m=field("Flowering period").apply(flowering)
values["Flowering onset month"]=m.str[0];values["Flowering duration months"]=m.str[1]
for k,src,pattern in [("Water dispersal strategy","Dispersal strategy","hydrochory"),("Evergreen leaves","Leaf life span","evergreen"),("Clonal herb","Growth form","clonal herb"),("Insect vector","Pollination syndrome","insect-pollination"),("Wind vector","Pollination syndrome","wind-pollination"),("Vegetative reproduction","Reproduction type","vegetatively")]:
    s=field(src).replace({"?":np.nan,"–":np.nan,"-":np.nan})
    values[k]=s.str.contains(pattern,case=False,regex=False).where(s.notna()).astype(float)
# Distinguish confirmed members of an explicit category from substring matches.
s=field("Growth form")
values["Clonal herb"]=s.str.split(r",\s*").apply(lambda z:float("clonal herb" in z) if isinstance(z,list) and z[0] not in ["–","-","?"] else np.nan)
keys=["site_region","island_band","site_band"]
pd.DataFrame(values).to_csv(cache/"compiled.csv.gz",index_label="sp")
rows=[]
for k,v in values.items():
    keep=w.assign(value=v).dropna(subset=["value"]).groupby(keys).candidate.transform("nunique")
    ids=keep.index[keep.eq(2)]
    row=compare(w,v,keys,np.random.default_rng(918))
    row.update(trait=k,known_reference=v[~w.candidate].notna().sum(),matched_reference=(~w.loc[ids,"candidate"]).sum())
    regions=sorted(w.loc[ids,"site_region"].unique())
    loo=[]
    for region in regions:
        sub=w.loc[ids].query("site_region!=@region")
        if sub.empty:continue
        dif=[];weights=[]
        for _,group in sub.groupby(keys):
            a=v[group.index[group.candidate]]; b=v[group.index[~group.candidate]]
            dif.append(a.mean()-b.mean());weights.append(len(a))
        loo.append(np.average(dif,weights=weights))
    row.update(regions=len(regions),loo_min=min(loo,default=np.nan),loo_max=max(loo,default=np.nan))
    rows.append(row)
res=pd.DataFrame(rows)
valid=res.p.notna()
res.loc[valid,"q"]=multipletests(res.loc[valid,"p"],method="fdr_bh")[1]
res.sort_values("p").to_csv(out/"traits.csv",index=False)
w[["species","candidate","url"]].assign(number_measured=pd.DataFrame(values).notna().sum(axis=1)).to_csv(out/"coverage.csv")
print(res.sort_values("p")[["trait","matched","matched_reference","difference","p","q","loo_min","loo_max"]].head(8).to_string(index=False))
