"""Compare list detection across decades in the same observer × Swiss grid × month."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import duckdb
import rasterio
from pyproj import Transformer
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

cache,lights,out=map(Path,sys.argv[1:4])
con=duckdb.connect(config={"threads":3,"memory_limit":"1600MB"})
con.execute(f"""create view records as select * exclude(species), case species when 'Coeloglossum viride' then 'Dactylorhiza viridis' when 'Anemonoides nemorosa' then 'Anemone nemorosa' else species end as species from read_parquet('{(cache/'eligible.parquet').as_posix()}')""")
con.sql("""select species, count(*) as eligible_records, min(year(day)) as first_year, max(year(day)) as last_year,
    sum((year(day)<=1999)::int) as early_records, sum((year(day)>=2020)::int) as recent_records
    from records where species in ('Orchis militaris','Dactylorhiza viridis','Lathyrus cicera') group by species"""
    ).df().to_csv(out/"focal_occurrence_coverage.csv",index=False)
stages=con.sql("""select lifeStage, count(*) as n from records group by lifeStage order by n desc""").df()
stages.to_csv(out/"field_stage_coverage.csv",index=False)
coords=con.sql("select distinct lat,lon from records").df()
forward=Transformer.from_crs(4326,2056,always_xy=True)
east,north=forward.transform(coords.lon,coords.lat)
coords["gx"]=(east//5000).astype(int)
coords["gy"]=(north//5000).astype(int)
con.register("coords",coords)
con.execute("""create temp table visits as
    select distinct day,gx,gy,lower(trim(recorder)) as observer,species
    from records join coords using(lat,lon)
    where (year(day) between 1992 and 1999 or year(day) between 2020 and 2024)
      and recorder is not null and length(trim(recorder))>=3
      and not regexp_matches(lower(recorder),'anonym|unknown|inconnu|unbekannt|not supplied')
      and not(month(day)=1 and dayofmonth(day)=1)""")
con.execute("""create temp table lists as select day,gx,gy,observer,month(day) as month,
    year(day) as year, (year(day)>=2020)::int as late,count(*) as richness
    from visits group by all having count(*)>=10""")
con.execute("""create temp table pairs as select gx,gy,observer,month from lists group by all
    having count(distinct late)=2""")
design=con.sql("""select * from lists join pairs using(gx,gy,observer,month)""").df()
if design.empty:raise ValueError("no repeated lists")
cells=design[["gx","gy"]].drop_duplicates().reset_index(drop=True)
xy=cells[["gx","gy"]].to_numpy()[:,None,:]*5000+np.array([(x,y) for x in range(500,5000,1000) for y in range(500,5000,1000)])
back=Transformer.from_crs(2056,4326,always_xy=True)
lon,lat=back.transform(xy[:,:,0].ravel(),xy[:,:,1].ravel())
irr=[]
for year in sorted(design.year.unique()):
    file=next(lights.glob(f"*{year}_{'calDMSP' if year<=2013 else 'simVIIRS'}.tif"))
    with rasterio.open(file) as src:
        r,c=map(np.asarray,rasterio.transform.rowcol(src.transform,lon,lat))
        window=rasterio.windows.Window.from_slices((r.min(),r.max()+1),(c.min(),c.max()+1))
        v=src.read(1,window=window)[r-r.min(),c-c.min()]
        assert v.min()>=0 and v.max()<=63
        irr.append(cells.assign(year=year,dn=v.reshape(-1,25).mean(axis=1)))
design=design.merge(pd.concat(irr),validate="m:1")
strata=["gx","gy","observer","month"]
design["sid"]=design.groupby(strata).ngroup()
con.register("design",design)
hits=con.sql("""select sid,late,species,count(*) as found from visits join design
    on visits.day=design.day and visits.gx=design.gx and visits.gy=design.gy and visits.observer=design.observer
    group by all""").df()
info=design.groupby(["sid","late"]).agg(lists=("day","size"),dn=("dn","mean"),richness=("richness","mean"),
                                         gx=("gx","first"),gy=("gy","first"))
a=info.xs(0,level="late").copy()
b=info.xs(1,level="late").reindex(a.index)
a["delta_dn"]=(b.dn-a.dn)/10
a["delta_effort"]=np.log(b.richness)-np.log(a.richness)
a["cell"]=a.gx*1000+a.gy
a["block"]=a.gx.floordiv(5)*1000+a.gy.floordiv(5)
cell=a.groupby("cell").agg(block=("block","first"),delta_dn=("delta_dn","mean"),initial=("dn","mean"),
                            change_richness=("delta_effort","mean"),pairs=("gx","size"))
rows=[]
for sp,g in hits.groupby("species"):
    row=dict(species=sp)
    z=(g.set_index(["sid","late"]).found/info.lists).dropna().reset_index(name="rate")
    z["signed"]=np.where(z.late.eq(1),1,-1)*z.rate*100
    total=z.groupby("sid").signed.sum().to_frame().join(a.cell).groupby("cell").signed.sum()
    s=cell.loc[total.index].copy()
    y=total/s.pairs
    for period,label in [(0,"early_pct"),(1,"recent_pct")]:
        counts=z.loc[z.late.eq(period),["sid","rate"]].set_index("sid").join(a.cell)
        per_cell=counts.groupby("cell").rate.sum().reindex(s.index,fill_value=0)
        row[label]=100*(per_cell/s.pairs).mean()
    assert np.isclose(row["recent_pct"]-row["early_pct"],y.mean())
    blocks=s.block.to_numpy()
    G=s.block.nunique()
    row.update(cells=len(s),blocks=G)
    rows.append(row)
    if G<10 or len(s)<20:continue
    diff=y.mean()
    scores=(y-diff).groupby(s.block).sum().to_numpy()
    se=np.sqrt((scores@scores)*G/(G-1))/len(s)
    if se<=0:continue
    row.update(delta_pp=diff,p_time=2*t.sf(abs(diff/se),G-1))
    matrix=s[["delta_dn","initial","change_richness"]].assign(y=y)
    matrix-=matrix.groupby(s.block).transform("mean")
    X=matrix.iloc[:,:3].to_numpy()
    Y=matrix.y.to_numpy()
    if np.linalg.matrix_rank(X)<3 or len(s)-G-3<=0:continue
    z=X[:,0]-X[:,1:]@np.linalg.lstsq(X[:,1:],X[:,0],rcond=None)[0]
    information=pd.Series(z*z).groupby(blocks).sum()
    count=(information>max(1e-12,information.sum()*1e-10)).sum()
    row["informative_light_blocks"]=count
    if count<10:continue
    inv=np.linalg.inv(X.T@X)
    beta=inv@X.T@Y
    score=pd.DataFrame(X*(Y-X@beta)[:,None]).groupby(blocks).sum().to_numpy()
    se=np.sqrt((inv@score.T@score@inv)[0,0]*G/(G-1)*(len(s)-1)/(len(s)-G-3))
    if se>0:
        row.update(pp_per_10dn=beta[0],p_light=2*t.sf(abs(beta[0]/se),G-1))
res=pd.DataFrame(rows)
for tag in ["time","light"]:
    p=res["p_"+tag].dropna()
    res.loc[p.index,"q_"+tag]=multipletests(p,method="fdr_bh")[1]
res.to_csv(out/"gbif_species.csv",index=False)
summary={"archive_rows":16302425,"eligible_occurrences":con.sql("select count(*) from records").fetchone()[0],
         "qualifying_lists":con.sql("select count(*) from lists").fetchone()[0],
         "paired_lists":len(design),"paired_strata":len(a),"paired_cells":len(cell),"paired_observers":design.observer.nunique(),
         "early_lists":int((design.late==0).sum()),"recent_lists":int((design.late==1).sum()),"blocks":int(cell.block.nunique())}
pd.Series(summary).rename("value").to_csv(out/"gbif_coverage.csv")
print(summary,flush=True)
print(res.query('q_time<.05 and delta_pp<0').nsmallest(12,'q_time').to_string(index=False),flush=True)
seen=con.sql("""select distinct year(day) as year,gx,gy,species from records join coords using(lat,lon)
  where year(day) between 1992 and 2014""").df()
seen.to_parquet(cache/"grid_species_year.parquet",index=False)




