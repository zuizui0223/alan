"""Evaluate absolute green-up to flowering separation in light-associated species."""
import sys
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests


def within(g,interaction):
    groups=[pd.factorize(pd.MultiIndex.from_frame(g[["source","unit_id"]]))[0],
            pd.factorize(pd.MultiIndex.from_frame(g[["source","year"]]))[0]]
    light=g.dn.to_numpy()/10
    x=np.column_stack((light,light*g.land_class.eq("island"))) if interaction else light[:,None]
    a=np.column_stack((x,g.gap))
    for _ in range(1000):
        maximum=0
        for codes in groups:
            counts=np.bincount(codes)
            for v in a.T:
                means=np.bincount(codes,weights=v)/counts
                maximum=max(maximum,max(abs(means)))
                v-=means[codes]
        if maximum<1e-8:break
    else:return {}
    x,y=a[:,:-1],a[:,-1]
    if np.linalg.matrix_rank(x,tol=1e-6)<x.shape[1]:return {}
    weight=np.linalg.pinv(x)
    beta=weight@y
    score=weight*(y-x@beta)
    blocks=pd.factorize(g.block)[0]
    n=len(y)
    df=n-sum(len(set(c)) for c in groups)-x.shape[1]+1
    if df<1:return {}
    out={}
    for i,label in enumerate(["light","difference"][:x.shape[1]]):
        info=np.bincount(blocks,weights=weight[i]**2)
        G=(info>1e-8*info.sum()).sum()
        if G<10:continue
        se=np.linalg.norm(np.bincount(blocks,weights=score[i]))*np.sqrt(G/(G-1)*(n-1)/df)
        if se<1e-10:continue
        ci=t.ppf(.975,G-1)*se
        out.update({label:beta[i],label+"_lo":beta[i]-ci,label+"_hi":beta[i]+ci,label+"_p":2*t.sf(abs(beta[i]/se),G-1)})
    return out


def main(cache,out):
    annual=pd.read_csv(out/"year_span20.csv.gz",low_memory=False)
    strict=pd.read_csv(out/"year_observed20.csv.gz",low_memory=False)
    trends=pd.read_csv(out/"long_term_screen.csv.gz",low_memory=False)
    monthly=pd.read_csv(out/"screen_span20.csv.gz",low_memory=False)
    strict_month=pd.read_csv(out/"screen_observed20.csv.gz",low_memory=False)
    primary=set(strict_month.loc[strict_month.q.lt(.05),"sp"])|set(monthly.loc[monthly.q.lt(.05),"sp"])|set(annual.loc[annual.q.lt(.05),"sp"])|set(strict.loc[strict.q.lt(.05),"sp"])|set(trends.loc[trends.q.lt(.05),"sp"])
    nominal=set(annual.loc[annual.p.lt(.05),"sp"]) if "p" in annual else set()
    if "p" in trends:nominal.update(trends.loc[trends.p.lt(.05),"sp"])
    targets=annual[annual.sp.isin(primary|nominal)|annual.species.eq("Plantago lanceolata")]
    field=pd.read_parquet(cache/"flowering.parquet")
    names=pd.read_csv(cache/"field_crosswalk.csv",dtype={"matched_key":str})
    lookup=annual.drop_duplicates("species",keep=False).set_index("species").sp
    names["sp"]=names.matched_key.where(names.matched_key.isin(set(annual.sp)),names.field_name.map(lookup))
    field=field.merge(names[["field_name","sp"]],left_on="species",right_on="field_name",validate="m:1")
    db=duckdb.connect(str(cache/"annual.duckdb"))
    db.register("ff",field)
    db.execute("""create or replace table field_observations as select f.*, exists (
      select 1 from year_surveys e join geo g using(site) join year_risk r using(site)
      where r.sp=f.sp and e.y=f.year and g.latitude=f.latitude and g.longitude=f.longitude
    ) as exact_survey_year from ff f""")
    field=db.sql("from field_observations").df()
    results=[]
    for row in targets.itertuples():
        f=field[field.sp.eq(row.sp)].copy()
        # Collapse aliases and duplicated exact station-years (one protocol per coordinate).
        f=f.sort_values("source").drop_duplicates(["unit_id","year"])
        r=dict(sp=row.sp,species=row.species,priority="FDR" if row.sp in primary else "exploratory",
               gbif_slope=getattr(row,"slope",np.nan),gbif_q=row.q,
               field_years=len(f),exact_survey_years=int(f.exact_survey_year.sum()))
        for land in ["island","mainland"]:
            side=f[f.land_class.eq(land)]
            r[land+"_coords"]=side.unit_id.nunique()
            r[land+"_landmasses"]=side.land_id.nunique()
            r[land+"_station_years"]=len(side)
            r[land+"_mean_gap_days"]=side.groupby("unit_id").gap.mean().mean()
            r[land+"_mean_signed_days"]=side.groupby("unit_id").signed_gap.mean().mean()
            repeated=side.groupby(["source","unit_id"]).filter(lambda s:s.year.nunique()>=3 and s.dn.nunique()>1) if len(side) else side
            r[land+"_variable_stations"]=repeated.unit_id.nunique()
            if repeated.unit_id.nunique()>=10 and repeated.block.nunique()>=10:
                for key,val in within(repeated,False).items():r[land+"_"+key]=val
        for group,key in [(["year"],"same_year"),(["source","year"],"same_source_year")]:
            if len(f):
                shared=f.groupby(group+["land_class"]).gap.mean().unstack()
                if {"island","mainland"}<=set(shared):
                    diff=(shared.island-shared.mainland).dropna()
                    r[key+"_comparisons"]=len(diff)
                    r[key+"_mean_difference_days"]=diff.mean()
        years=f.groupby("year").land_class.nunique()
        matched=f[f.year.isin(years[years.eq(2)].index)]
        if len(matched):matched=matched.groupby(["source","unit_id"]).filter(lambda s:len(s)>=3 and s.dn.nunique()>1)
        r["common_years"]=matched.year.nunique()
        rep=matched.groupby("land_class").agg(coords=("unit_id","nunique"),blocks=("block","nunique"))
        if (len(rep)==2 and rep.coords.min()>=10 and rep.blocks.min()>=10 and
            matched.loc[matched.land_class.eq("island"),"land_id"].nunique()>=5):
            r.update({k:v for k,v in within(matched,True).items() if k.startswith("difference")})
        results.append(r)
    d=pd.DataFrame(results)
    for col in ["island_light_p","mainland_light_p","difference_p"]:
        if col not in d:d[col]=np.nan
        ok=d[col].notna()
        d[col+"_q"]=np.nan
        if ok.any():d.loc[ok,col+"_q"]=multipletests(d.loc[ok,col],method="fdr_bh")[1]
    d.to_csv(out/"candidate_flowering.csv",index=False)
    names.to_csv(out/"flowering_name_crosswalk.csv",index=False)
    print("Field response candidates",len(d),"island and mainland",sum(d.island_coords.gt(0)&d.mainland_coords.gt(0)),flush=True)


if __name__=="__main__":main(*map(Path,sys.argv[1:3]))
