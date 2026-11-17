"""Effect modification across the same wet-plant survey sites."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from alan_wet_shore import frozen, out as shore

out=shore.parent/"wet_drivers"
out.mkdir(exist_ok=True)
g=frozen("coordinate_trends.csv.gz").merge(pd.read_csv(shore/"distance.csv").dropna(subset="distance_m"),on="unit_id",validate="1:1")
g["near"]=1-g.distance_m/1000
g["light10"]=g.light/10
controls=["light10","light0","near","area0","freq0","area_trend","freq_trend","temp0","rain0","temp_trend","rain_trend","built0","built","effort","days","spring","summer","autumn","start","span","years","island","nc","nr"]
modifiers=["near","area0","area_trend","freq_trend","temp0","rain0","temp_trend","rain_trend","built0","built","light0","island"]
def fit(data,modifier):
    f=data.groupby("block").filter(lambda z:len(z)>1 and np.ptp(z.light)>1e-8).copy()
    scale=1. if modifier in ["near","island"] else f[modifier].quantile(.75)-f[modifier].quantile(.25)
    # Always test the interaction even when the empirical IQR vanishes.
    if scale<=1e-9:scale=f[modifier].std()
    if scale<=1e-9:return {}
    f["interaction"]=f.light10*(f[modifier]-f[modifier].median())/scale
    names=["interaction",*controls,"rate_difference"]
    arr=(f[names]-f.groupby("block")[names].transform("mean")).to_numpy()
    x=arr[:,:-1];y=arr[:,-1]
    rank=np.linalg.matrix_rank(x)
    if np.linalg.matrix_rank(x[:,1:])==rank:return {}
    inv=np.linalg.pinv(x)
    coef=inv@y
    v=inv[0]
    score=pd.Series(v*(y-x@coef)).groupby(f.block.values).sum()
    share=pd.Series(v*v).groupby(f.block.values).sum();share/=share.sum()
    n=len(f); blocks=len(score);k=(share>1e-9).sum()
    se=np.linalg.norm(score)*np.sqrt(k/(k-1)*(n-1)/(n-blocks-rank))
    half=t.ppf(.975,k-1)*se
    return dict(sites=n,blocks=blocks,scale=scale,estimate=coef[0],low=coef[0]-half,high=coef[0]+half,p=2*t.sf(abs(coef[0]/se),k-1),effective_blocks=1/(share@share))
rows=[]
for label,data in [("all",g),("brightening",g[g.light>1e-8])]:
    for m in modifiers:
        rows.append(dict(scope=label,modifier=m,**fit(data,m)))
res=pd.DataFrame(rows)
for _,s in res.dropna(subset=["p"]).groupby("scope"):
    res.loc[s.index,"q"]=multipletests(s.p,method="fdr_bh")[1]
res.to_csv(out/"site_interactions.csv",index=False)
print(res[res.scope.eq("all")].sort_values("p").round(3).to_string(index=False))
