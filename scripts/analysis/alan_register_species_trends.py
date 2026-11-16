"""Separate within-cell temporal trends and year-adjusted light associations in recorded shares."""
from pathlib import Path
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
import sys
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

cache, wcvp, pr85, ledger, out = map(Path,sys.argv[1:6])
sys.path.insert(0,str(pr85.resolve()))
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits

sites=pd.read_parquet(cache/"long_sites.parquet")
m=pd.read_parquet(cache/"long_members.parquet")
sites["cell"]=sites.y*7201+sites.x
sites["block"]=sites.y.floordiv(100)*73+sites.x.floordiv(100)
m=m.merge(sites[["year","x","y","cell"]].assign(idx=sites.index),validate="many_to_one")
site_rows={c:g.index.to_numpy() for c,g in sites.groupby("cell")}
year=sites.year.to_numpy()
effort=sites.effort.to_numpy()
light=sites.light.to_numpy()
cell=sites.cell.to_numpy()
block=sites.block.to_numpy()

def estimate(outcome, predictor, controls, cell_codes, block_codes):
    # Each sampled cell has unit total weight. Absorb fixed site effects exactly.
    counts=np.bincount(cell_codes)
    weight=1/counts[cell_codes]
    X=np.column_stack([outcome,predictor,controls])
    for col in X.T:
        col-= (np.bincount(cell_codes,weights=col)/counts)[cell_codes]
    X*=np.sqrt(weight[:,None])
    Y=X[:,:2].copy()
    C=X[:,2:]
    Y-=C@np.linalg.lstsq(C,Y,rcond=1e-10)[0]
    y,z=Y.T
    ss=z@z
    if ss<1e-7:return np.nan,np.nan
    b=(z@y)/ss
    G=block_codes.max()+1
    scores=np.bincount(block_codes,weights=z*(y-b*z))
    df=len(z)-len(counts)-np.linalg.matrix_rank(C)-1
    if df<=0:return np.nan,np.nan
    se=np.sqrt(scores@scores*G/(G-1)*(len(z)-1)/df)/ss
    if se<=0:return b,np.nan
    return b,2*t.sf(abs(b/se),G-1)

rows=[]
for sp,g in m.groupby("sp"):
    cc=g.cell.unique()
    nb=len(np.unique(block[g.idx]))
    row=dict(sp=sp,positive_cells=len(cc),blocks=nb,observation_uuids=g.n.sum())
    rows.append(row)
    if len(cc)<10 or nb<10:continue
    idx=np.concatenate([site_rows[c] for c in cc])
    ccode=pd.factorize(cell[idx])[0]
    bcode=pd.factorize(block[idx])[0]
    y=pd.Series(g.n.to_numpy(),index=g.idx).reindex(idx,fill_value=0).to_numpy()*100/effort[idx]
    years=year[idx]
    nuisance=np.log(effort[idx])
    row["time_pp_per_decade"],row["p_time"]=estimate(y,(years-2000)/10,nuisance,ccode,bcode)
    ycat=np.eye(years.max()-years.min()+1)[years-years.min(),1:]
    row["light_pp_per_10dn"],row["p_light"]=estimate(y,light[idx]/10,np.column_stack([nuisance,ycat]),ccode,bcode)
d=pd.DataFrame(rows).set_index("sp")
for tag in ["time","light"]:
    p=d["p_"+tag].dropna()
    d.loc[p.index,"q_"+tag]=multipletests(p,method="fdr_bh")[1]
# Intersection-union test for joint directional declines; adjust across all jointly estimable species.
p_time=d.p_time.div(2).where(d.time_pp_per_decade.lt(0),1-d.p_time.div(2))
p_light=d.p_light.div(2).where(d.light_pp_per_10dn.lt(0),1-d.p_light.div(2))
p=pd.concat([p_time,p_light],axis=1).dropna().max(axis=1)
d.loc[p.index,"q_joint_negative"]=multipletests(p,method="fdr_bh")[1]
names=pd.read_csv(next(wcvp.glob("*accepted_native_distribution*")),usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"]).drop_duplicates()
assert names.wcvp_accepted_plant_name_id.is_unique
d["species"]=names.set_index("wcvp_accepted_plant_name_id").wcvp_accepted_name.reindex(d.index)
# Fill only if the WCVP resolution ledger labels the source binomial itself accepted.
lookup=pd.read_csv(next(wcvp.glob("*resolution*")))
extra=lookup[lookup.wcvp_name_status_seen.eq("accepted") & lookup.wcvp_match_state.eq("unique_exact_name_to_accepted_id")]
extra=extra[["wcvp_accepted_plant_name_id","source_species_exact_binomial"]].drop_duplicates()
extra=extra[~extra.wcvp_accepted_plant_name_id.duplicated(keep=False)].set_index("wcvp_accepted_plant_name_id")
d["species"]=d.species.fillna(extra.source_species_exact_binomial.reindex(d.index))
d=d.join(direct_plant_traits(ledger),on="species")
traits=pd.read_csv(out.parent/"joint_night_autonomy/species.csv").set_index("species")
d=d.join(traits[["night_share"]],on="species")
share=pd.read_csv(out.parent/"species_interpretation/island_sign_reproductive_evidence.csv").set_index("plant")
d=d.join(share[["transport_review_status","glopl_median_log_response_ratio","island_pp","mainland_pp"]],on="species")
d.to_csv(out/"species.csv")
for tag,flag in {"estimated":d.q_time.notna(),"time_down":d.q_time.lt(.05)&d.time_pp_per_decade.lt(0),"joint":d.q_joint_negative.lt(.05)}.items():
    s=d[flag]
    print(tag,len(s),s.autonomous_selfing_capacity.value_counts().to_dict(),"night data",s.night_share.count())
print(d.nsmallest(10,"q_joint_negative")[["species","time_pp_per_decade","light_pp_per_10dn","q_joint_negative"]])
