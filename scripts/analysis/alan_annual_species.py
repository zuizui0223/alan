"""Screen both directions in annual site-month records, retaining surveyed non-detections."""
from pathlib import Path
import sys
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests


def fit(g):
    n=len(g)
    unit=g.unit_id.nunique()
    row=dict(cells=n,coordinates=unit,years=g.y.nunique(),blocks=g.block.nunique(),
             island_coordinates=g.loc[g.land_class.eq("island"),"unit_id"].nunique(),
             source_sites=g.site.nunique(),survey_days=int(g.visits.sum()),detection_days=int(g.hits.sum()),
             zero_month_years=int(g.hits.eq(0).sum()),post_last_detection_site_years=int(
                 len(g.loc[g.y>g.last_seen,["site","y"]].drop_duplicates())))
    if unit<20 or row["blocks"]<10:return dict(row,status="insufficient_spatial_replication")
    strata=pd.factorize(pd.MultiIndex.from_frame(g[["site","m"]]))[0]
    years=pd.factorize(g.y)[0]
    # Each site-month-year receives equal weight. Remove FE by alternating projections.
    rates=g.hits.to_numpy()/g.visits.to_numpy()
    fractions=[c for c in g if c.startswith("month_fraction_")]
    a=np.column_stack((g.dn/10,np.log(g.richness),g[fractions].to_numpy(),rates*100))
    for iteration in range(2000):
        change=0
        for codes in (strata,years):
            counts=np.bincount(codes)
            for v in a.T:
                means=np.bincount(codes,weights=v)/counts
                change=max(change,max(abs(means)))
                v-=means[codes]
        if change<1e-8:break
    else:return dict(row,status="fixed_effects_not_converged")
    light,controls,y=a[:,0],a[:,1:-1],a[:,-1]
    x=light-controls@np.linalg.lstsq(controls,light,rcond=1e-9)[0]
    ss=x@x
    if ss<1e-8:return dict(row,status="no_adjusted_light_contrast")
    blocks=pd.factorize(g.block)[0]
    info=np.bincount(blocks,weights=x*x)
    informative=int((info>max(1e-10,ss*1e-9)).sum())
    row.update(informative_blocks=informative,largest_block_light_share=info.max()/ss)
    if informative<10:return dict(row,status="insufficient_informative_blocks")
    b=x@y/ss
    u=y-b*light
    u-=controls@np.linalg.lstsq(controls,u,rcond=1e-9)[0]
    df=n-len(np.unique(strata))-len(np.unique(years))-np.linalg.matrix_rank(controls)
    if df<1:return dict(row,status="no_residual_df")
    scores=np.bincount(blocks,weights=x*u)/ss
    se=np.linalg.norm(scores)*np.sqrt(informative/(informative-1)*(n-1)/df)
    if se<1e-10:return dict(row,status="degenerate_variance")
    ci=t.ppf(.975,informative-1)*se
    return dict(row,status="estimated",slope=b,lo=b-ci,hi=b+ci,p=2*t.sf(abs(b/se),informative-1))


def main(cache,out,annual=False):
    db=duckdb.connect(str(cache/"annual.duckdb"),config={"threads":3,"memory_limit":"1600MB"})
    # Include all surveyed years of each source-site and month in which a species
    # was ever detected. Ordering is symmetric: zeros may precede or follow records.
    db.execute("""create or replace table supported as select * from eligible semi join
       (select site,m from eligible group by all having count(distinct y)>=2) using(site,m)""")
    db.execute("""create or replace table risk as select distinct h.sp,h.site,h.m
       from hits h semi join supported s using(site,y,m)""")
    db.execute("""create or replace table last_seen as select site,sp,max(y) last_seen from hits group by all""")
    roster=db.sql("""select n.sp,n.species,count(distinct h.site) as survey_sites
       from species_names n left join hits h using(sp) group by all order by n.sp""").df()
    years=db.sql("select site,count(distinct y) as sample_years from eligible group by site").df()
    db.register("sy",years)
    if annual:
        parts=",".join(f"sum(if(m={m},visits,0))/sum(visits) as month_fraction_{m}" for m in range(1,12))
        db.execute("""create or replace table year_surveys as select site,y,unit_id,block,land_id,land_class,dn,
            sum(visits) as visits,sum(richness*visits)/sum(visits) as richness,"""+parts+
            " from eligible group by all")
        db.execute("""create or replace table year_hits as select h.site,h.y,h.sp,sum(h.hits) as hits
             from hits h semi join year_surveys using(site,y) group by all""")
        db.execute("create or replace table year_risk as select distinct site,sp from year_hits")
    rows=[]
    for bucket in range(48):
        if not annual:g=db.sql(f"""select r.sp,s.*,coalesce(h.hits,0) as hits,sy.sample_years,l.last_seen
          from risk r join supported s using(site,m) join sy using(site) join last_seen l using(site,sp)
          left join hits h using(site,y,m,sp) where hash(r.sp)%48={bucket}""").df()
        if annual:
            g=db.sql(f"""select r.sp,s.*,0 as m,coalesce(h.hits,0) as hits,sy.sample_years,l.last_seen
              from year_risk r join year_surveys s using(site) join sy using(site) join last_seen l using(site,sp)
              left join year_hits h using(site,y,sp) where hash(r.sp)%48={bucket}""").df()
        for sp,grp in g.groupby("sp",sort=False):
            for cohort,sub in [("span20",grp),("observed20",grp[grp.sample_years>=20])]:
                if len(sub):rows.append(dict(sp=sp,cohort=cohort,**fit(sub)))
        if bucket%8==0:print("species partition",bucket,len(rows),flush=True)
    estimates=pd.DataFrame(rows)
    output=[]
    for cohort in ["span20","observed20"]:
        d=roster.merge(estimates[estimates.cohort.eq(cohort)].drop(columns="cohort"),how="left",validate="1:1")
        d["status"]=d.status.fillna("no_eligible_long_survey_station" if annual else "no_comparable_months_or_sites")
        d["q"]=np.nan
        if "p" in d:
            ix=d.p.notna()
            if ix.any():d.loc[ix,"q"]=multipletests(d.loc[ix,"p"],method="fdr_bh")[1]
        d["direction"]="unresolved"
        if "slope" in d:
            for flag,label in [(d.slope.lt(0),"negative"),(d.slope.gt(0),"positive")]:
                d.loc[flag & d.q.lt(.05),"direction"]=label
        if annual:d=d.rename(columns={"zero_month_years":"surveyed_zero_years"})
        prefix="year" if annual else "screen"
        d.to_csv(out/f"{prefix}_{cohort}.csv.gz",index=False,compression="gzip")
        output.append(d[d.q.lt(.05)].assign(cohort=cohort))
        print(cohort,d.status.value_counts().to_dict(),d.direction.value_counts().to_dict(),flush=True)
    pd.concat(output).to_csv(out/("annual_candidates.csv" if annual else "light_associated_species.csv"),index=False)
    db.execute("checkpoint")

if __name__=="__main__":
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main(*map(Path,sys.argv[1:3]),annual="--all-visits" in sys.argv[3:])
