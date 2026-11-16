"""Species-specific changes at 20-year repeat coordinates; select joint decline/brightening evidence."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import duckdb
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

def cluster_mean(y,blocks):
    n=len(y)
    codes,G=pd.factorize(blocks)[0],len(np.unique(blocks))
    if G<10 or n<20:return None
    b=np.mean(y)
    scores=np.bincount(codes,weights=y-b)
    se=np.linalg.norm(scores)*np.sqrt(G/(G-1))/n
    if se<1e-10:return None
    return b,se,G

def model(s):
    if len(s)<20 or s.block.nunique()<10:return None
    z=s[["delta","light0","rich_change","change"]].copy()
    z["delta"]/=10
    mat=(z-z.groupby(s.block).transform("mean")).to_numpy()
    y=mat[:,-1]
    controls=mat[:,1:3]
    q=mat[:,0]-controls@np.linalg.lstsq(controls,mat[:,0],rcond=None)[0]
    ss=q@q
    if ss<1e-8:return None
    codes=pd.factorize(s.block)[0]
    info=np.bincount(codes,weights=q*q)
    informative=np.count_nonzero(info>max(1e-10,ss*1e-9))
    if informative<10:return None
    b=(q@y)/ss
    error=y-b*mat[:,0]
    error-=controls@np.linalg.lstsq(controls,error,rcond=None)[0]
    scores=np.bincount(codes,weights=q*error)/ss
    df=len(s)-len(info)-1-np.linalg.matrix_rank(controls)
    if df<1:return None
    adj=np.sqrt(informative/(informative-1)*(len(s)-1)/df)
    scores*=adj
    se=np.linalg.norm(scores)
    if se<1e-10:return None
    return b,se,informative,dict(zip(pd.unique(s.block),scores))

def output_fit(row,tag,fit):
    if fit is None:return
    b,se,G=fit[:3]
    crit=t.ppf(.975,G-1)
    row.update({tag:b,tag+"_lo":b-crit*se,tag+"_hi":b+crit*se,tag+"_pneg":t.cdf(b/se,G-1)})

def adjusted(d,src,dst):
    d[dst]=np.nan
    if src not in d:return
    p=d[src].dropna()
    if len(p):d.loc[p.index,dst]=multipletests(p,method="fdr_bh")[1]

def run(cache,out):
    db=duckdb.connect(str(cache/"global_gbif.duckdb"),config={"threads":3,"memory_limit":"1400MB"})
    sql="""select p.sp,g.unit_id as xy,g.land_class,g.land_id,g.block,avg(g.delta) as delta,avg(g.light0) as light0,
     avg(g.rich_change) as rich_change,avg(p.recent-p.early) as change from pairs p join site_exposure g using(site)
     where g.land_class in ('island','mainland') group by all"""
    db.execute("create or replace table plant_coordinate_change as "+sql)
    names=db.sql("select specieskey, species from names").df().set_index("specieskey").species
    rows=[]
    for bucket in range(16):
        table=db.sql(f"select * from plant_coordinate_change where hash(sp)%16={bucket}").df()
        for sp,g in table.groupby("sp"):
            bright=g[g.delta.gt(0)]
            row=dict(sp=sp,species=names[sp],coordinates=len(g),island_coordinates=g.land_class.eq("island").sum(),
                     mainland_coordinates=g.land_class.eq("mainland").sum(),blocks=g.block.nunique(),brightened=len(bright),
                     mean_change=g.change.mean(),brightness_change_iqr=g.delta.quantile(.75)-g.delta.quantile(.25))
            if len(g)>=20 and g.block.nunique()>=10:
                output_fit(row,"slope",model(g))
                output_fit(row,"bright_change",cluster_mean(bright.change.to_numpy(),bright.block.to_numpy()))
            rows.append(row)
    d=pd.DataFrame(rows)
    for col in ["slope_pneg","bright_change_pneg"]:
        if col not in d:d[col]=np.nan
    p=d[["slope_pneg","bright_change_pneg"]]
    d["p_joint"]=p.max(axis=1,skipna=False)
    adjusted(d,"p_joint","q_joint");adjusted(d,"slope_pneg","q_slope")
    d.to_csv(out/"species_screen.csv",index=False)
    hits=d[d.q_joint.lt(.05)]
    extra=[]
    exploratory=d[d.p_joint.notna() & d.slope.lt(0) & d.bright_change.lt(0)] if d.p_joint.notna().any() else d.iloc[:0]
    for r in exploratory.itertuples():
        g=db.execute("select * from plant_coordinate_change where sp=?",[r.sp]).df()
        a=g[g.land_class.eq("island")];b=g[g.land_class.eq("mainland")]
        row=dict(sp=r.sp,species=r.species,global_joint_q=r.q_joint,island_coords=len(a),mainland_coords=len(b),islands=a.land_id.nunique(),
                 island_blocks=a.block.nunique(),mainland_blocks=b.block.nunique())
        left=model(a) if a.land_id.nunique()>=5 else None
        right=model(b)
        output_fit(row,"island",left);output_fit(row,"mainland",right)
        if left and right:
            G=min(left[2],right[2])
            scores=pd.Series(left[3]).subtract(pd.Series(right[3]),fill_value=0)
            se=np.linalg.norm(scores)
            diff=left[0]-right[0]
            row.update(difference=diff,lo=diff-t.ppf(.975,G-1)*se,hi=diff+t.ppf(.975,G-1)*se,
                       p=2*t.sf(abs(diff/se),G-1))
        extra.append(row)
    dif=pd.DataFrame(extra,columns=list(dict.fromkeys(["sp","species","island_coords","mainland_coords","islands",
         "island_blocks","mainland_blocks","island","mainland","difference","lo","hi","p"]+[k for r in extra for k in r])))
    adjusted(dif,"p","q")
    dif.to_csv(out/"selected_land_contrasts.csv",index=False)
    print("species, gradient eligible, joint eligible, selected",len(d),d.slope_pneg.count(),d.p_joint.count(),len(hits),flush=True)
    if len(hits):print(hits.nsmallest(10,"q_joint")[["species","slope","bright_change","q_joint"]].to_string(index=False))

if __name__=="__main__":run(*map(Path,sys.argv[1:3]))
