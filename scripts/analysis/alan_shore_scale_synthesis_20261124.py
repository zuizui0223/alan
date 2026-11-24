"""Replay placement/spectrum contrasts and water-distance threshold sensitivity.
Run from ALAN repo root, with scipy/pandas/numpy and the fixed blobs available.
Python adaptation of executed in-session JS; no original workbook re-extraction.
"""
import gzip
import io
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t

BLOBS = {
 "units":"c544c4af5c39677f8d5d4bacbd40ef2a518868fb",
 "taxa":"d99b0e97d6608d8420342394ca5e9c65a89011c5",
 "coords":"1c025bcf6c28a6a3cd125ad05628169e126aef0d",
 "distance":"c2bd22ab52c4674fedfb39a12f67cf847a5e2c3d"
}
def read_blob(key):
    b = subprocess.check_output(["git","show",BLOBS[key]])
    if b[:2] == b"\x1f\x8b":
        b = gzip.decompress(b)
    return pd.read_csv(io.BytesIO(b))

def placement():
    u, tax = read_blob("units"), read_blob("taxa")
    assert u.unit_id.is_unique and len(u)==30
    group = tax.groupby("unit_id")
    aqua = tax[tax.origin_as_published_figure3.eq("aquatic")].groupby("unit_id")["count"].sum()
    terra = tax[tax.origin_as_published_figure3.eq("terrestrial")].groupby("unit_id")["count"].sum()
    u = u.set_index("unit_id")
    u["aquatic"] = aqua.reindex(u.index,fill_value=0)
    u["terrestrial"] = terra.reindex(u.index,fill_value=0)
    u["all"] = u.aquatic + u.terrestrial
    assert (u.aquatic==u.aquatic_origin_capture).all()
    assert (u.terrestrial==u.terrestrial_origin_capture).all()
    assert u["all"].sum()==4699
    def taxa(names):
        return tax[tax.family_source_spelling.isin(names)].groupby("unit_id")["count"].sum().reindex(u.index,fill_value=0)
    for k in ["caenidae","leptoceridae","chironomidae"]:
        u[k]=taxa([k])
    u["aquatic_ex_caen"]=u.aquatic-u.caenidae
    u["aquatic_ex_two"]=u.aquatic-u.caenidae-u.leptoceridae
    u["terrestrial_ex_eph"]=u.terrestrial-taxa(["ephydridae"])
    u["trichoptera"]=taxa(["leptoceridae","hydroptillidae","rhyacophilidae","psychomiidae","polycentropodidae","phryganeidae","helicopsyche","hydropsychidae"])
    colours=["green","purple","red","white"]
    menus={"all":colours,"no_violet":["green","red","white"],"green_white":["green","white"],**{c:[c] for c in colours}}
    endpoints=["all","aquatic","terrestrial","aquatic_ex_caen","aquatic_ex_two","caenidae","leptoceridae","chironomidae","trichoptera"]
    terms={e:{e:1} for e in endpoints}
    terms["origin_interaction"]={"aquatic":1,"terrestrial":-1}
    terms["origin_interaction_ex_eph"]={"aquatic":1,"terrestrial_ex_eph":-1}
    rows=[]
    for name,menu in menus.items():
        for endpoint,term in terms.items():
            hrr={}; vv=[]; zero=False
            for h in ["littoral","riparian"]:
                for arm,cs in [("lit",menu),("dark",["dark"])]:
                    cells=[u[(u.habitat==h)&(u.colour==c)] for c in cs]
                    assert all(len(z)==3 for z in cells)
                    mu=np.mean([z[list(term)].mean().values for z in cells],axis=0)
                    if (mu<=0).any():
                        zero=True;break
                    hrr[h,arm]=sum(s*np.log(m) for s,m in zip(term.values(),mu))
                    for z in cells:
                        score=(z[list(term)].values @ (np.array(list(term.values()))/mu))/len(cs)
                        vv.append(np.var(score,ddof=1)/len(z))
                if zero:break
            if zero:
                rows.append({"menu":name,"endpoint":endpoint,"status":"zero_mean"});continue
            a=hrr["littoral","lit"]-hrr["littoral","dark"]
            b=hrr["riparian","lit"]-hrr["riparian","dark"]
            theta=a-b; v=sum(vv); df=v*v/sum(x*x/2 for x in vv) if v>0 else np.inf
            hw=t.ppf(.975,df)*np.sqrt(v)
            rows.append({"menu":name,"endpoint":endpoint,"status":"estimated","ratio":np.exp(theta),"variance":v,"df":df,"lower":np.exp(theta-hw),"upper":np.exp(theta+hw),"littoral_RR":np.exp(a),"riparian_RR":np.exp(b)})
    return pd.DataFrame(rows)

def macro():
    g=read_blob("coords").merge(read_blob("distance").dropna(subset=["distance_m"]),on="unit_id",validate="1:1")
    assert len(g)==1706
    env=["light0","temp0","rain0","temp_trend","rain_trend","built0","built","effort","days","spring","summer","autumn","start","span","years","island","nc","nr","area0","freq0","area_trend","freq_trend"]
    outcomes=["candidate_rate","reference_rate","rate_difference"]
    g["light10"]=g.light/10
    g["linear_near"]=1-g.distance_m/1000
    g["linear_near_light"]=g.linear_near*g.light10
    result=[]
    for scope,frame in [("all",g),("brightened",g[g.light>1e-8])]:
        f=frame.groupby("block").filter(lambda z:len(z)>1 and z.light.max()-z.light.min()>1e-8).copy()
        n=len(f); kblock=f.block.nunique()
        for cut,extra in [("linear",False)]+[(c,e) for e in (False,True) for c in (40,60,100)]:
            f["near"]=f.linear_near if cut=="linear" else (f.distance_m<=cut).astype(float)
            f["near_light"]=f.near*f.light10
            pred=["near_light","light10","near",*env]+(["linear_near","linear_near_light"] if extra else [])
            dm=f[pred+outcomes]-f.groupby("block")[pred+outcomes].transform("mean")
            X=dm[pred].to_numpy(); Y=dm[outcomes].to_numpy()
            rank=np.linalg.matrix_rank(X); inv=np.linalg.pinv(X); df=n-kblock-rank
            a=np.zeros(len(pred));a[0]=1
            if not np.allclose(a,a@inv@X):continue
            v=a@inv; beta=inv@Y; b=a@beta; resid=Y-X@beta
            scores=pd.DataFrame(v[:,None]*resid).groupby(f.block.values).sum()
            shares=pd.Series(v*v).groupby(f.block.values).sum();shares/=shares.sum()
            k=int((shares>1e-9).sum())
            se=np.sqrt((scores**2).sum().values*k/(k-1)*(n-1)/df)
            half=t.ppf(.975,k-1)*se
            for i,o in enumerate(outcomes):
                result.append(dict(scope=scope,cut=cut,add_continuous=extra,outcome=o,estimate=b[i],low=b[i]-half[i],high=b[i]+half[i],sites=n,blocks=kblock,effective_blocks=1/(shares@shares),near_sites=None if cut=="linear" else int(f.near.sum())))
    out=pd.DataFrame(result)
    for scope,target in [("all",3.8953320748238482),("brightened",1.8884812937990745)]:
        val=out[(out.scope==scope)&(out.cut=="linear")&(out.outcome=="rate_difference")].estimate.iloc[0]
        assert abs(val-target)<1e-7
    return out
if __name__=="__main__":
    out=Path("results/all_species_forest/biotic_interactions/shore_scale_replay")
    out.mkdir(exist_ok=True,parents=True)
    placement().to_csv(out/"placement.csv",index=False)
    macro().to_csv(out/"macro.csv",index=False)
