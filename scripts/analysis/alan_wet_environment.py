"""Secondary analysis of the frozen species cohort; one row per plant species."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

repo = Path(__file__).resolve().parents[2]
root = repo / "results/all_species_forest"
out = root / "wet_environment"
out.mkdir(exist_ok=True)
cols = ["light", "temperature", "precipitation", "built"]
raw = pd.read_parquet(repo.parent/"gbif_annual/all_species_rf_sites.parquet", columns=["sp","unit_id",*cols])
assert not raw.duplicated(["sp","unit_id"]).any()
f = pd.read_csv(root/"decline_traits/moisture_evidence.csv.gz").join(raw.groupby("sp")[cols].median(),on="sp")
assert len(f)==3901 and f.candidate.sum()==283
keys=["site_region","island_band","site_band"]
f=f.dropna(subset=["indicator",*keys,*cols]).copy()
g=f.groupby(keys).candidate
f=f[g.transform("min").ne(g.transform("max"))].copy()
g=f.groupby(keys).candidate
n=g.transform("sum")
f["weight"]=np.where(f.candidate,1,n/(g.transform("size")-n))
f["candidate"]=f.candidate.astype(int)
region=[]
for k,s in f.groupby("site_region"):
    c=s[s.candidate.eq(1)]
    b=s[s.candidate.eq(0)]
    region.append(dict(region=k,candidates=len(c),references=len(b),difference=c.indicator.mean()-np.average(b.indicator,weights=b.weight)))
regional=pd.DataFrame(region).sort_values("candidates",ascending=False)
regional.to_csv(out/"regions.csv",index=False)
# Quantiles and scaling are fixed at the full matched sample for delete-region fits.
scales=[]
for col in cols:
    q1,med,q3=f[col].quantile([.25,.5,.75])
    scales.append(dict(term=col,q25=q1,median=med,q75=q3))
    f[col]=(f[col]-med)/(q3-q1)
    f["by_"+col]=f[col]*f.candidate
def fit(s,terms):
    mat=s[["indicator",*terms]].copy()
    by=[s[k] for k in keys]
    means=mat.mul(s.weight,axis=0).groupby(by).transform("sum").div(s.weight.groupby(by).transform("sum"),axis=0)
    x=(mat-means).mul(np.sqrt(s.weight),axis=0)
    coeff,_,rank,_=np.linalg.lstsq(x[terms],x.indicator,rcond=None)
    assert rank==len(terms)
    return coeff
rows=[]
for model,terms in [("frozen",["candidate"]),("environment",["candidate",*cols]),
                    ("moderation",["candidate",*cols,*["by_"+col for col in cols]])]:
    beta=fit(f,terms)
    leave=np.array([fit(f[f.site_region.ne(k)],terms) for k in regional.region])
    size=len(leave)
    half=t.ppf(.975,size-1)*np.sqrt((size-1)/size*((leave-leave.mean(axis=0))**2).sum(axis=0))
    for term,est,h in zip(terms,beta,half):
        rows.append(dict(model=model,term=term,estimate=est,low=est-h,high=est+h,
                         candidates=int(f.candidate.sum()),references=int((1-f.candidate).sum()),regions=size))
effects=pd.DataFrame(rows)
ix=effects.term.str.startswith("by_")
se=(effects.loc[ix,"high"]-effects.loc[ix,"low"])/(2*t.ppf(.975,len(regional)-1))
pvals=2*t.sf(abs(effects.loc[ix,"estimate"]/se),len(regional)-1)
effects.loc[ix,"holm_p_four_interactions"]=multipletests(pvals,method="holm")[1]
effects.to_csv(out/"environment.csv",index=False)
pd.DataFrame(scales).to_csv(out/"exposure_scale.csv",index=False)
assert np.isclose(rows[0]["estimate"],0.5922104455502862)
