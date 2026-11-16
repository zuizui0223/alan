"""Pooled equal-species interval estimates, with optional GBIF grid-year confirmation."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import t
from alan_swiss_flower_gap import summarize

def pooled(panel):
    panel=panel.groupby(["station","species"]).filter(lambda g:len(g)>=3).copy()
    panel=panel.groupby("species").filter(lambda g:g.station.nunique()>=10 and g.block.nunique()>=10)
    outcomes=["abs_gap","signed_gap","flower","sos"]
    v=panel[["dn",*outcomes]].to_numpy(float)
    v[:,0]/=10
    w=1/panel.groupby("species").station.transform("size").to_numpy(float)
    within=[pd.factorize(pd.MultiIndex.from_frame(panel[keys]))[0] for keys in [["species","station"],["species","year"]]]
    for _ in range(500):
        before=v.copy()
        for codes in within:
            denom=np.bincount(codes,weights=w)
            for col in v.T:
                col-=(np.bincount(codes,weights=w*col)/denom)[codes]
        if np.max(np.abs(v-before))<1e-9:break
    assert np.max(np.abs(v-before))<1e-8
    z=v[:,0]
    G,cl=np.unique(panel.block,return_inverse=True)
    n=len(G)
    ss=(w*z)@z
    # Each species panel is connected across years; use conservative absorbed df.
    df=len(v)-sum(int(g.max())+1 for g in within)
    rows=[]
    for j,tag in enumerate(outcomes,1):
        y=v[:,j]
        b=(w*z)@y/ss
        s=np.bincount(cl,weights=w*z*(y-b*z))
        se=np.sqrt(s@s/ss**2*n/(n-1)*(len(v)-1)/df)
        ci=se*t.ppf(.975,n-1)
        rows.append(dict(endpoint=tag,rows=len(v),species=panel.species.nunique(),stations=panel.station.nunique(),
                         years=panel.year.nunique(),blocks=n,estimate=b,lo=b-ci,hi=b+ci))
    return pd.DataFrame(rows)

if __name__=="__main__":
    cache,out=map(Path,sys.argv[1:3])
    gbif=cache/"grid_species_year.parquet"
    results=[]
    coverage=[]
    for phase,file in [("onset","aligned.parquet"),("50pct","full_flower_aligned.parquet")]:
        d=pd.read_parquet(cache/"pheno"/file)
        coverage.append(d.groupby("species").agg(stations=("station","nunique"),rows=("year","size"),first=("year","min"),last=("year","max")).assign(phase=phase,cohort="all"))
        results.append(pooled(d).assign(phase=phase,cohort="all"))
        if gbif.exists():
            seen=pd.read_parquet(gbif,filters=[("species","in",list(d.species.unique()))])
            d=d.merge(seen,on=["year","gx","gy","species"],validate="m:1")
            coverage.append(d.groupby("species").agg(stations=("station","nunique"),rows=("year","size"),first=("year","min"),last=("year","max")).assign(phase=phase,cohort="gbif_grid_year"))
            summarize(d).to_csv(out/f"gbif_overlap_{phase}.csv",index=False)
            results.append(pooled(d).assign(phase=phase,cohort="gbif_grid_year"))
    pd.concat(coverage).to_csv(out/"phase_coverage.csv")
    r=pd.concat(results)
    r.to_csv(out/"pooled_intervals.csv",index=False)
    print(r[r.endpoint.eq("abs_gap")].to_string(index=False))

